#!/usr/bin/env python
from gi.repository import Ufo

#{{{ parameters
scale            = 2
number_of_images = 1

ring_start     = 5 / scale
ring_end       = 7 / scale
ring_step      = 2  / scale
ring_count     = ( ring_end - ring_start )  / ring_step + 1
ring_thickness = 4 / scale
#img_path = '/home/chuan/DATA/upiv/Image0.tif'; xshift = 0; yshfit = 200
#img_path = '/home/chuan/DATA/upiv/60_Luft_frame50.tif'; xshift = 150; yshift = 0
#img_path = 'data/image1.tif'; xshift=0; yshift=0
img_path = 'data/input/SampleC'; xshift=0; yshift=0
res_tif = './ordfilt.tif'
res_txt = './res.txt'
#}}}
#{{{ filters
filters = {
    'read': {
       'property': { 
           'path': img_path,
           'number': number_of_images
       } 
    },
    'crop': {
        'property': {
            'x': xshift, 'y': yshift, 
            'width': 1024, 'height':1024  
        }
    },
    'duplicater': {
        'name': 'buffer',
        'property': {'dup_count': ring_count}
    },
    'denoise': {
        'property': {'matrix_size': 26/scale}
    },
    'contrast': {
        'property': {'remove_high': 0}
    },
    'ordfilt': {},
    'ring_pattern': {
        'name': 'of-ring_pattern',
        'property': {
            'method': 1,
            'start': ring_start,
            'end': ring_end,
            'step': ring_step,
            'thickness': ring_thickness,
            'width': 1024/scale,
            'height': 1024/scale
        }
    },
    'ring_loop': {
        'name': 'buffer',
        'property': {'loop': 1}
    },
    'filter_particle': {
        'name': 'filter-particle'
    },
    'ring_writer': {},
    'write': {
        'property': { 'filename': res_tif }
    },
    'reduce': {},
    'stack': {
        'property': {'number': ring_count}
    }
}
#}}}


pm = Ufo.PluginManager()
for k,v in filters.iteritems():
    name = v.get('name') or k
    prop = v.get('property') or {}
    locals()[k] = pm.get_task(name)
    locals()[k].set_properties(**prop)


g = Ufo.TaskGraph()
# read image and pre-processing
g.connect_nodes(read, crop)
if scale == 2:
    g.connect_nodes(crop, reduce)
    g.connect_nodes(reduce, denoise)
else:
    g.connect_nodes(crop, denoise)
g.connect_nodes(denoise, contrast)
# Ordered Search Filter
g.connect_nodes(ring_pattern, ring_loop)
g.connect_nodes_full(contrast, ordfilt, 0)
g.connect_nodes_full(ring_loop, ordfilt, 1)
g.connect_nodes(ordfilt, filter_particle)
g.connect_nodes(filter_particle, ring_writer)


import sys; sys.exit()

sched = Ufo.Scheduler()
sched.set_properties(enable_tracing=False)
sched.run(g)
