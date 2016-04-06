import argparse
import os
import math
from subprocess import call

#### setup

parser = argparse.ArgumentParser()
parser.add_argument('--challenge', required=True, choices=['1A', '1B', '1C', '2A', '2B', '3A', '3B'])
parser.add_argument('--data', choices=['tiny', 'half', 'full', '1k', '5k', '10k', '20k', '50k', '100k'])
parser.add_argument('--approx', nargs=2, type=float) # 2 args (fraction, iterations)
args = parser.parse_args()

if args.challenge[0] != '1' and args.data == None:
	print 'INPUT ERROR: you forgot to specify --data'
	quit()

print args

if args.approx != None:
	approx_fraction = args.approx[0]
	approx_iterations = args.approx[1]
	if approx_fraction <= 0.0 or approx_fraction >= 1.0:
		print 'INPUT ERROR: first argument to --approx is the sub-sampling fraction; this value [x] should be 0.0 < x < 1.0'
		quit()
	approx_iterations = int(math.floor(approx_iterations))
	print approx_fraction, approx_iterations

#### work

SCRIPT_PATH = '../smc_het_eval/SMCScoring.py'
DATA_ROOT_DIR = '../../data/kbuckets'
OUTPUT_DIR = './output-%s' % args.challenge[0]

try:
	os.mkdir(OUTPUT_DIR)
except OSError as e:
	if 'File exists' in e:
		pass
	else:
		print e

command_pieces = []

## python call

command_pieces.append('python')
command_pieces.append(SCRIPT_PATH)

## challenge

command_pieces.append('-c')
command_pieces.append(args.challenge)

## pred/truth/vcf

if args.challenge[0] == '1':
	command_pieces.append('--predfiles')
	command_pieces.append('%s/%s/pred.txt' % (DATA_ROOT_DIR, args.challenge))
	command_pieces.append('--truthfiles')
	command_pieces.append('%s/%s/truth.txt' % (DATA_ROOT_DIR, args.challenge))

	if args.challenge == '1C':
		command_pieces.append('--vcf')
		command_pieces.append('%s/%s/scoring.vcf' % (DATA_ROOT_DIR, args.challenge))

else:
	filetype = 'txt' if args.challenge[1] == 'A' else 'txt.gz'

	if args.challenge[0] == '2':
		command_pieces.append('--predfiles')
		command_pieces.append('%s/%s/%s/pred.%s' % (DATA_ROOT_DIR, args.challenge, args.data, filetype))
		command_pieces.append('--truthfiles')
		command_pieces.append('%s/%s/%s/truth.%s' % (DATA_ROOT_DIR, args.challenge, args.data, filetype))
	elif args.challenge[0] == '3':
		command_pieces.append('--predfiles')
		command_pieces.append('%s/%s/%s/pred.2%s.%s' % (DATA_ROOT_DIR, args.challenge, args.data, args.challenge[1],filetype))
		command_pieces.append('%s/%s/%s/pred.3%s.%s' % (DATA_ROOT_DIR, args.challenge, args.data, args.challenge[1],filetype))
		command_pieces.append('--truthfiles')
		command_pieces.append('%s/%s/%s/truth.2%s.%s' % (DATA_ROOT_DIR, args.challenge, args.data, args.challenge[1],filetype))
		command_pieces.append('%s/%s/%s/truth.3%s.%s' % (DATA_ROOT_DIR, args.challenge, args.data, args.challenge[1],filetype))

	command_pieces.append('--vcf')
	command_pieces.append('%s/%s/%s/scoring.vcf' % (DATA_ROOT_DIR, args.challenge, args.data))

## output

command_pieces.append('-o')
command_pieces.append('%s/%s_score.txt' % (OUTPUT_DIR, args.challenge))

## approx

if args.approx != None:
	command_pieces.append('--approx')
	command_pieces.append(str(approx_fraction))
	command_pieces.append(str(approx_iterations))

## done

command = ' '.join(command_pieces)
print command

## the call

call(command, shell=True)