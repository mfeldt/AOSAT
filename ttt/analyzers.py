
import aosat
from aosat.analyze import frg_analyzer, sps_analyzer, run, setup
import os
example_path = os.path.join(os.path.split(aosat.__file__)[0],'examples')
example_file = os.path.join(example_path,'example_analyze_closed_loop.setup')
aosat.aosat_cfg.CFG_SETTINGS = aosat.aosat_cfg.configure(example_file)
aosat.aosat_cfg.CFG_SETTINGS["totalnumber"]=1
sd=setup()
analyzers=[sps_analyzer(sd)]
run(analyzers)
