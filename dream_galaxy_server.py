#!/usr/bin/env python

import os
import argparse
import tornado
import tornado.web
import requests
import thread
import threading
import time
import json
import subprocess
from jinja2 import Template

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORK_DIR = "/tmp/dream"

class GalaxyProxy(tornado.web.RequestHandler):
    
    def galaxy_init(self, galaxy, apikey, keyfile):
        if keyfile is not None:
            with open(keyfile) as handle:
                apikey = handle.read().rstrip()
        self.remote = RemoteGalaxy(galaxy, apikey)
        self.galaxy = galaxy
        self.apikey = apikey
        


class MainHandler(GalaxyProxy):
    
    def initialize(self, galaxy, apikey, keyfile, **kwds):
        with open(os.path.join(TEMPLATE_DIR, "main.html")) as handle:
            self.page = Template(handle.read())
        
        self.galaxy_init(galaxy, apikey, keyfile)

    def get(self):
        workflows = self.remote.workflow_list()
        
        host_url = "%s://%s" % (self.request.protocol, self.request.host)
        login=False
        if 'login' not in self.request.arguments:
            login=True
        self.write(self.page.render(host_url=host_url, workflows=workflows, login=login))

class FormHandler(GalaxyProxy):

    def initialize(self, galaxy, apikey, submitter, keyfile, **kwds):
        self.submitter = submitter
        with open(os.path.join(TEMPLATE_DIR, "info_form.html")) as handle:
            self.page = Template(handle.read())
        self.galaxy_init(galaxy, apikey, keyfile)


    def get(self):
        if 'workflow' not in self.request.arguments:
            with open(os.path.join(TEMPLATE_DIR, "submit.html")) as handle:
                page = Template(handle.read())
            self.write(page.render(message="<h1>Missing Workflow Selection</h1>", dest="."))
        else:
            workflow = self.request.arguments['workflow'][0]
            self.write(self.page.render(workflow=workflow))
                

class ValidateHandler(GalaxyProxy):

    def initialize(self, galaxy, apikey, submitter, keyfile, **kwds):
        self.submitter = submitter
        with open(os.path.join(TEMPLATE_DIR, "submit.html")) as handle:
            self.page = Template(handle.read())
        self.galaxy_init(galaxy, apikey, keyfile)

    def get(self):
        workflow_path = self.request.arguments['workflow'][0]
        workflow = self.remote.get(workflow_path + "/download")
        if not self.submitter.running:
            meta = {}
            for k,v in self.request.arguments.items():
                meta[k] = v[0]
            self.submitter.submission = {
                'workflow' : self.galaxy + workflow_path + "/download",
                'apikey' : self.apikey,
                'flags' : "--no-upload --check", 
                'meta' : {},
                'ok_message' : "Workflow Looks OK"
            }
            message = "<h1>Checking</h1> %s" % (json.dumps(self.submitter.submission))
        else:
            message = "Already Working on submission"
        self.write(self.page.render(message=message, dest="monitor"))
        
class SubmitHandler(GalaxyProxy):

    def initialize(self, galaxy, apikey, submitter, keyfile, **kwds):
        self.submitter = submitter
        with open(os.path.join(TEMPLATE_DIR, "submit.html")) as handle:
            self.page = Template(handle.read())
        self.galaxy_init(galaxy, apikey, keyfile)

    def post(self):
        workflow_path = self.request.arguments['workflow'][0]
        workflow = self.remote.get(workflow_path + "/download")
        if not self.submitter.running:
            meta = {}
            for k,v in self.request.arguments.items():
                meta[k] = v[0]
            self.submitter.submission = {
                'workflow' : self.galaxy + workflow_path + "/download",
                'apikey' : self.apikey,
                'synapse_email' : self.request.arguments['synapse_email'][0],
                'synapse_apikey' : self.request.arguments['synapse_apikey'][0],   
                'meta' : meta,
                'flags' : "--no-upload",
                'ok_message' : "Workflow Submitted"
            }
            message = "<h1>Submitting</h1> %s" % (json.dumps(self.submitter.submission))
        else:
            message = "Already Working on submission"
        self.write(self.page.render(message=message, dest="monitor"))
                
class MonitorHandler(tornado.web.RequestHandler):

    def initialize(self, submitter, **kwds):
        self.submitter = submitter
        with open(os.path.join(TEMPLATE_DIR, "monitor.html")) as handle:
            self.page = Template(handle.read())

    def get(self):
        self.write(self.page.render(log=self.submitter.log))


class Submitter(threading.Thread):
    
    def run(self):
        self.quit = False
        self.log = ""
        self.running = False
        self.submission = None
        while not self.quit:
            time.sleep(1)
            if self.submission is not None:
                self.running = True
                
                if not os.path.exists(WORK_DIR):
                    os.mkdir(WORK_DIR)
                
                with open(os.path.join(WORK_DIR, "submission.json"), "w") as handle:
                    handle.write(json.dumps(self.submission['meta']))
                
                cmd_line_template = "{submit_cmd} \
--meta {meta_path} \
--synapse_email {synapse_email} \
--synapse_key {synapse_key} \
--meta {meta_path} \
--apikey {apikey} \
--workdir {workdir} \
--workflow {workflow} \
{flags}" 
                
                cmd_line = cmd_line_template.format(
                    submit_cmd=os.path.join(BASE_DIR, "dream_galaxy_submit"),
                    meta_path=os.path.join(WORK_DIR, "submission.json"),
                    synapse_email=self.submission.get('synapse_email', "test@test.com"),
                    synapse_key=self.submission.get('synapse_apikey', "NA"),
                    apikey=self.submission['apikey'],
                    workdir=WORK_DIR,
                    flags=self.submission['flags'],
                    workflow=self.submission['workflow']
                )
                
                self.log = "Running: %s\n" % (cmd_line)
                proc = subprocess.Popen(cmd_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
                def watch_and_log(stream):
                    while True:
                        line = stream.readline()
                        if not line:
                            break
                        self.log += line
                thread.start_new_thread(watch_and_log, (proc.stdout,))
                thread.start_new_thread(watch_and_log, (proc.stderr,))
                proc.wait()
                if proc.returncode == 0:
                    self.log += "\n%s" % (self.submission['ok_message'])
                else:
                    self.log += "\nFailed"                    
                self.running = False
                self.submission = None


class RemoteGalaxy(object):

    def __init__(self, url, api_key, path_mapping={}):
        self.url = url
        self.api_key = api_key
        self.path_mapping = path_mapping

    def get(self, path, params = {}):
        c_url = self.url + path
        params['key'] = self.api_key
        req = requests.get(c_url, params=params)
        return req.json()

    def post(self, path, payload, params={}):
        c_url = self.url + path
        params['key'] = self.api_key
        logging.debug("POSTING: %s %s" % (c_url, json.dumps(payload)))
        req = requests.post(c_url, data=json.dumps(payload), params=params, headers = {'Content-Type': 'application/json'} )
        print req.text
        return req.json()

    def post_text(self, path, payload, params=None):
        c_url = self.url + path
        if params is None:
            params = {}
        params['key'] = self.api_key
        logging.debug("POSTING: %s %s" % (c_url, json.dumps(payload)))
        req = requests.post(c_url, data=json.dumps(payload), params=params, headers = {'Content-Type': 'application/json'} )
        return req.text

    def download_handle(self, path):
        url = self.url + path
        logging.info("Downloading: %s" % (url))
        params = {}
        params['key'] = self.api_key
        r = requests.get(url, params=params, stream=True)
        return r

    def download(self, path, dst):
        r = self.download_handle(path)
        dsize = 0L
        if hasattr(dst, 'write'):
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    dst.write(chunk)
                    dsize += len(chunk)
        else:
            with open(dst, "wb") as handle:
                for chunk in r.iter_content(chunk_size=1024):
                    if chunk:
                        handle.write(chunk)
                        handle.flush()
                        dsize += len(chunk)
        logging.info("Downloaded: %s bytes" % (dsize))

    def create_library(self, name):
        lib_create_data = {'name' : name}
        library = self.post('/api/libraries', lib_create_data)
        library_id = library['id']
        return library_id

    def library_find(self, name):
        for d in self.library_list():
            if d['name'] == name:
                return d
        return None

    def library_list(self):
        return self.get("/api/libraries")

    def library_list_contents(self, library_id):
        return self.get("/api/libraries/%s/contents" % library_id)

    def library_find_contents(self, library_id, name):
        for a in self.library_list_contents(library_id):
            if a['name'] == name:
                return a
        return None

    def library_get_contents(self, library_id, ldda_id):
        return self.get("/api/libraries/%s/contents/%s" % (library_id, ldda_id))

    def get_hda(self, history, hda):
        return self.get("/api/histories/%s/contents/%s" % (history, hda))

    def get_dataset(self, id, src='hda' ):
        return self.get("/api/datasets/%s?hda_ldda=%s" % (id, src))

    def download_hda(self, history, hda, dst):
        meta = self.get_hda(history, hda)
        self.download(meta['download_url'], dst)

    def history_list(self):
        return self.get("/api/histories")

    def get_history(self, history):
        return self.get("/api/histories/%s" % (history))

    def get_history_contents(self, history):
        return self.get("/api/histories/%s/contents?details=all" % (history))

    def get_provenance(self, history, hda, follow=False):
        if follow:
            return self.get("/api/histories/%s/contents/%s/provenance" % (history, hda), {"follow" : True})
        else:
            return self.get("/api/histories/%s/contents/%s/provenance" % (history, hda))

    def add_workflow(self, wf):
        self.post("/api/workflows/upload", { 'workflow' : wf } )

    def workflow_list(self):
        return self.get("/api/workflows")

    def get_workflow(self, wid):
        return self.get("/api/workflows/%s" % (wid))

    def call_workflow(self, request):
        return self.post("/api/workflows", request, params={'step_details' : True} )

    def get_job(self, jid):
        return self.get("/api/jobs/%s" % (jid), {'full' : True} )

    def library_paste_file(self, library_id, library_folder_id, name, datapath, uuid=None, metadata=None):
        datapath = os.path.abspath(datapath)
        found = False
        for ppath, dpath in self.path_mapping.items():
            if datapath.startswith(ppath):
                datapath = os.path.join(dpath, os.path.relpath(datapath, ppath))
                found = True
                break
        if not found:
            raise Exception("Path not in mounted lib_data directories: %s" % (datapath))
        data = {}
        data['folder_id'] = library_folder_id
        data['file_type'] = 'auto'
        data['name'] = name
        if uuid is not None:
            data['uuid'] = uuid
        data['dbkey'] = ''
        data['upload_option'] = 'upload_paths'
        data['create_type'] = 'file'
        data['link_data_only'] = 'link_to_files'
        if metadata is not None:
            data['extended_metadata'] = metadata
        data['filesystem_paths'] = datapath
        logging.info("Pasting %s: %s" % (name, datapath))
        libset = self.post("/api/libraries/%s/contents" % library_id, data)
        print libset
        return libset[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy", default="http://localhost")
    parser.add_argument("--apikey", default=None)
    parser.add_argument("--keyfile", default=None)
    parser.add_argument("--port", type=int, default=8080)

    args = parser.parse_args()

    submitter = Submitter()
    submitter.start()
    
    config = {
        'galaxy' : args.galaxy, 
        'apikey' : args.apikey, 
        'keyfile' : args.keyfile, 
        'submitter' : submitter
    }
    
    application = tornado.web.Application([
        (r"/", MainHandler, config),
        (r"/info_form", FormHandler, config),
        (r"/submit", SubmitHandler, config),
        (r"/validate", ValidateHandler, config),
        (r"/monitor", MonitorHandler, config),

    ])

    application.listen(args.port)
    tornado.ioloop.IOLoop.instance().start()


