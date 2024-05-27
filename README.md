# bounce_inflation

the chains can be downloaded at this URL https://pan.cstcloud.cn/s/IZcjxeWQtk. The Chinese character “下载” means "Download".

The montepython package can be found at https://github.com/brinckmann/montepython_public.

The CLASS package can be found at https://github.com/lesgourg/class_public. One should replace the original file class_public_version/python/class.pyx by the class.pyx file in this repository and recompile the CLASS package.

The SPT-3G data we used in this work can be downloaded from https://github.com/ksardase/SPT3G-montepython.  A new version can be found at https://github.com/SouthPoleTelescope/spt3g_y1_dist.

One can compile the primordial power spectrum file of bounce-inflation by command:

`g++ spectrum.cpp -o power_spectrum_bounce_inflation -lpthread`
