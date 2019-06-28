# libresample16
16 bits resample library inspired by https://ccrma.stanford.edu/~jos/resample/Free_Resampling_Software.html

- API interface separated from filter build and sndfile
- sample interleave instead of using 2 separated buffers
- limited memory copy (create small memory buffer between calls)
- buildfilter produces .h directly
