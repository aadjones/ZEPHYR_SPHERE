LDFLAGS_COMMON = -framework Accelerate -framework GLUT -framework OpenGL -L/opt/local/lib/ -L/usr/local/lib/x86_64 -L/usr/local/lib -lstdc++ -lpng -lz -ljpeg -msse2 -lfftw3 -lm
CFLAGS_COMMON = -v -c -Wall -I../../src/integrators -I../../src/linearalgebra -I../../src/geometry -I../../src/util -I../../src/glvu -DDO_PNG_OUT=0 -I./ -I/usr/local/include -I../../ -I../../src/Eigen/ -O3 -DNO_FFT -fopenmp -msse2 -lstdc++ -I/opt/local/include/

#LDFLAGS_COMMON = -framework Accelerate -framework GLUT -framework OpenGL ${GLVUFLAGS} -L/usr/local/lib/x86_64 -L/usr/local/lib -lstdc++ -lpng -lz -ljpeg -fopenmp -msse2 -fprofile-arcs -ftest-coverage -lfftw3 -lm
#CFLAGS_COMMON = -c -Wall -I../../src/integrators -I../../src/linearalgebra -I../../src/geometry -I../../src/util -I../../src/glvu -DDO_PNG_OUT=0 -I${GLVU_INCLUDE} -I./ -I/usr/local/include -I../../ -I../../src/Eigen/ -O1 -DNO_FFT -fopenmp -msse2 -fprofile-arcs -ftest-coverage
