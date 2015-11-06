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



class MainHandler(tornado.web.RequestHandler):
    
    def initialize(self, galaxy, apikey):
        with open(os.path.join(TEMPLATE_DIR, "main.html")) as handle:
            self.page = Template(handle.read())
        self.remote = RemoteGalaxy(galaxy, apikey)

    def get(self):
        workflows = self.remote.workflow_list()
        
        host_url = "%s://%s" % (self.request.protocol, self.request.host)
        login=False
        if 'login' not in self.request.arguments:
            login=True
        self.write(self.page.render(host_url=host_url, workflows=workflows, login=login))

class FormHandler(tornado.web.RequestHandler):

    def initialize(self, galaxy, apikey, submitter):
        self.submitter = submitter
        self.galaxy = galaxy
        self.apikey = apikey
        self.remote = RemoteGalaxy(galaxy, apikey)
        with open(os.path.join(TEMPLATE_DIR, "info_form.html")) as handle:
            self.page = Template(handle.read())

    def get(self):
        workflow = self.request.arguments['workflow'][0]
        self.write(self.page.render(workflow=workflow))
                

class SubmitHandler(tornado.web.RequestHandler):

    def initialize(self, galaxy, apikey, submitter):
        self.submitter = submitter
        self.galaxy = galaxy
        self.apikey = apikey
        self.remote = RemoteGalaxy(galaxy, apikey)
        with open(os.path.join(TEMPLATE_DIR, "submit.html")) as handle:
            self.page = Template(handle.read())

    def post(self):
        workflow_path = self.request.arguments['workflow'][0]
        workflow = self.remote.get(workflow_path + "/download")
        if not self.submitter.running:
            self.submitter.submission = {
                'workflow' : self.galaxy + workflow_path + "/download",
                'apikey' : self.apikey,
                'synapse_email' : self.request.arguments['synapse_email'][0],
                'synapse_apikey' : self.request.arguments['synapse_apikey'][0],                
            }
        self.write(self.page.render(message=json.dumps(self.request.arguments)))
                
class MonitorHandler(tornado.web.RequestHandler):

    def initialize(self, galaxy, apikey, submitter):
        self.submitter = submitter
        self.galaxy = galaxy
        self.apikey = apikey
        self.remote = RemoteGalaxy(galaxy, apikey)
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
                
                cmd_line_template = "./dream_galaxy_submit \
--user {user} \
--password {password} \
--team-name {team_name} \
--name {name} \
--apikey {apikey} \
--workflow {workflow} --no-upload" 
                
                cmd_line = cmd_line_template.format(
                    user="test",
                    name="test_workflow",
                    password="password",
                    team_name="team-name",
                    apikey=self.submission['apikey'],
                    workflow=self.submission['workflow']
                )
                
                self.log = "Running: %s" % (cmd_line)
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
                    self.log += "\nSubmission Done"
                else:
                    self.log += "\nSubmission Failed"                    
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
    parser.add_argument("--apikey", default="admin")    
    parser.add_argument("--port", type=int, default=8080)

    args = parser.parse_args()

    submitter = Submitter()
    submitter.start()

    application = tornado.web.Application([
        (r"/", MainHandler, {'galaxy' : args.galaxy, 'apikey' : args.apikey}),
        (r"/info_form", FormHandler, {'galaxy' : args.galaxy, 'apikey' : args.apikey, 'submitter' : submitter}),
        (r"/submit", SubmitHandler, {'galaxy' : args.galaxy, 'apikey' : args.apikey, 'submitter' : submitter}),
        (r"/monitor", MonitorHandler, {'galaxy' : args.galaxy, 'apikey' : args.apikey, 'submitter' : submitter}),

    ])

    application.listen(args.port)
    tornado.ioloop.IOLoop.instance().start()


