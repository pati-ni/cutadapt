#!/usr/bin/env python
import cProfile as profile
import pstats
import sys

p = pstats.Stats(sys.argv[1])
p.sort_stats('time').print_stats(30)
#from cutadapt.scripts.cutadapt import main
#profile.run('main()', 'cutadapt.prof')

#p = pstats.Stats('cutadapt.prof')
#p.sort_stats('cumulative').print_stats(20)
#p.sort_stats('time').print_stats(50)
