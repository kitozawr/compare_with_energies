# compare_with_energies

Compare data obtained at one computer with pulse energies stored at another computer. Computer timers are not sinchronized with each other. 
Takes into account skips in data recording.
The method relies on correlation between the pulse energy and signal intensity.

## Dependencies
Use functions from data_proc module: https://github.com/d-push/filamentation_data_proc

Requires python3 to be installed on local machine. Depends on numpy, matplolib, scipy, os, sys, shutil, struct and argparse python modules.
