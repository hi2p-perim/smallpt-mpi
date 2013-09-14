# smallpt-mpi

Hybrid MPI / OpenMP parallel implementation of
[smallpt](http://www.kevinbeason.com/smallpt/).

## Comparison

### smallpt (original version)

Configuration:

* OS : Debian 7
* Processor : Core i7-3870X

Commands and results:
> $ g++ -O3 -fopenmp -o smallpt smallpt.cpp  
> $ time ./smallpt 1000  
> Rendering (1000 spp) 100.00%  
> real    **3m0.072s**  
> user    35m39.082s  
> sys     0m0.032s

### smallpt-mpi

Experiment on Beowulf cluster
with 1-2 gen. Intel Core i7 processors
and gigabit ethernet interconnection.

* OS : Debian 7
* MPI implementation : OpemMPI
* Processors
  - Master node
     - Celeron G1610
  - Computation nodes
     - Core i7-880
     - Core i7-2600
     - Core i7-3930K
     - Core i7-3970X
     - Xeon e3-1230 v2

Commands and results:
> $ mkdir build  
> $ cd build  
> $ cmake -DCMAKE\_BUILD\_TYPE=Release ..  
> $ make  
> $ time mpirun -bynode -np 6 ./smallpt-mpi 1000  
> real    **1m12.508s**  
> user    0m6.996s  
> sys     1m4.976s

## License

The MIT License (MIT)

Copyright (c) 2013 Hisanari Otsu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

