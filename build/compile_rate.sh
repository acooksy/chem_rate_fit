#!/bin/bash
rate_root_dir=/home/acooksy/rate/
cd $rate_root_dir
gfortran src/rate-main.f src/rate-res.f src/uncert.f src/studentt.f include/xmrqmin.f include/gaussj.f -o rate

