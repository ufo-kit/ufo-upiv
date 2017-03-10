#!/usr/bin/env python
# useage:
#   ./piv-hough.py [config filename]
# default config file is config.py
import sys
import os
import logging
import argparse
import ConfigParser
from upiv.pivjob import PivJob

logfile = '.hough.log'
logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode='w', 
                    format='%(name)s %(levelname)s %(message)s')

def main():
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('configuration', nargs='?', default='config_hough.cfg')
    parser.add_argument('-g', '--graph')
    parser.add_argument('-p', '--profiling', action='store_true', default=False)
    args = parser.parse_args(sys.argv[1:])

    # parse configuration file
    cf_parser = ConfigParser.ConfigParser()
    if cf_parser.read(args.configuration) != []:
        parms = {k: eval(v) for k,v in cf_parser.items('parms')}
        config = {k: eval(v) for k,v in cf_parser.items('config')}
    else:
        print "Unable to open configuration file", args.configuration
        sys.exit(-1)

    # overwrite config/parms with command line arguments
    if args.graph:
        config['graph'] = args.graph
    if args.profiling:
        config['profiling'] = args.profiling

    # run pivjob
    j = PivJob(parms)
    j.profiling = config.get('profiling') or False
    j.deviceCPU = config.get('deviceCPU') or False
    j.schedfixed = config.get('schedfixed') or True
    j.run(config.get('graph') or 'azimu')

class PivJob(PivJob):
    def setup_tasks(self):
        self.setup_basic_tasks()
        p = self.parms
        sc = self.parms.scale

        self.add_task('crop', x=p.xshift, y=p.yshift, width=p.width, height=p.height)
        if 1:
            self.add_task('contrast', 'piv_contrast', 
                    c1=p.contrast_c1, c2=p.contrast_c2, gamma=p.contrast_gamma)
        else:
            self.add_task('contrast', remove_high=p.remove_high)
        self.add_task('rescale', factor=1.0/sc)
        self.add_task('input_fft', 'fft', dimensions=2)
        self.add_copy_task('bc_image') 

        self.add_task('ring_fft', 'fft', dimensions=2)
        self.add_task('ring_stack', 'stack', number=p.ring_number)
        self.add_task('ring_loop', 'loop', count=p.number)
        self.add_task('ring_convolution', 'complex_mult')
        self.add_task('ring_slice', 'slice')
        self.add_task('ring_pattern', 
                start=p.ring_start/sc, end=p.ring_end/sc, step=p.ring_step/sc, 
                thickness=p.ring_thickness/sc, method=p.ring_method, 
                width=p.width/sc, height=p.height/sc)

        self.add_task('ifft', dimensions=2)
        self.add_task('likelihood', 'hough-likelihood', 
                masksize=p.likelihoodmask, maskinnersize=p.likelihoodmaskinner,
                threshold=p.thld_likely)
        self.add_task('likelihood-image', 'hough-likelihood-image', 
                masksize=p.likelihoodmask, maskinnersize=p.likelihoodmaskinner)
        self.add_task('cand', 'candidate-filter', 
                ring_start=p.ring_start, ring_step=p.ring_step,
                ring_end=p.ring_end, scale=p.scale )
        self.add_task('azimu', 'azimuthal-test', thread=p.threads,
                azimu_thld = p.thld_azimu, likelihood_thld = p.thld_likely)

    def setup_graph(self, flag):
        if flag==0 or flag=='azimu':
            b1 = self.branch('read', 'crop', 'contrast', 'bc_image', 'rescale', 'input_fft')
            b2 = self.branch('ring_pattern', 'ring_stack', 'ring_fft', 'ring_loop')
            b3 = self.branch('bc_image')
            b4 = self.branch('ring_convolution', 'ifft', 'likelihood', 'cand')
            b5 = self.branch('azimu', 'ring_writer')
            self.graph.merge_branch(b1, b2, b4)
            self.graph.merge_branch(b3, b4, b5)
        elif flag==1 or flag=='candidate':
            b1 = self.branch('read', 'crop', 'contrast', 'rescale', 'input_fft')
            b2 = self.branch('ring_pattern', 'ring_stack', 'ring_fft', 'ring_loop')
            b3 = self.branch('ring_convolution', 'ifft', 'likelihood', 'cand', 'ring_writer')
            self.graph.merge_branch(b1, b2, b3)
        elif flag==2 or flag=='likelihood':
            b1 = self.branch('read', 'crop', 'contrast', 'rescale', 'input_fft')
            b2 = self.branch('ring_pattern', 'ring_stack', 'ring_fft', 'ring_loop')
            b3 = self.branch('ring_convolution', 'ifft', 'likelihood-image', 'write')
            self.graph.merge_branch(b1, b2, b3)
        elif flag==3 or flag=='hough':
            b1 = self.branch('read', 'crop', 'contrast', 'rescale', 'input_fft')
            b2 = self.branch('ring_pattern', 'ring_stack', 'ring_fft', 'ring_loop')
            b3 = self.branch('ring_convolution', 'ifft', 'write')
            self.graph.merge_branch(b1, b2, b3)
        elif flag==4 or flag=='contrast':
            b1 = self.branch('read', 'crop', 'contrast', 'write')
        elif flag==5 or flag=='azimu2':
            b1 = self.branch('read', 'crop', 'bc_image', 'contrast', 'rescale', 'input_fft')
            b2 = self.branch('ring_pattern', 'ring_stack', 'ring_fft', 'ring_loop')
            b3 = self.branch('bc_image')
            b4 = self.branch('ring_convolution', 'ifft', 'likelihood', 'cand')
            b5 = self.branch('azimu', 'ring_writer')
            self.graph.merge_branch(b1, b2, b4)
            self.graph.merge_branch(b3, b4, b5)
        else:
            print flag
            pass

if __name__ == '__main__':
    main()
