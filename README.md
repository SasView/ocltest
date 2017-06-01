# Minimal OpenCL program on Windows

This is a sample project to demonstrate a minimal OpenCL project on Windows. To compile it you only need 
MinGW64 and a graphics driver with an OpenCL runtime. Then compile it with

	gcc -I. main.c C:\Windows\System32\OpenCL.dll -o main.exe

No SDK or anything else needed. More details in [this article](http://arkanis.de/weblog/2014-11-25-minimal-opencl-development-on-windows).

Run using:

    main.exe


With tinycc x86 on windows use:

    tcc.exe -I. main.c c:/Windows/system32/OpenCL.dll

TinyCC is available from https://github.com/sasview/tinycc.git.  If this is installed as a sister repository, then use:

    ../tinycc/tinycc/x86/tcc.exe -I../tinycc/tinycc/include -I. main.c c:/Windows/system32/OpenCL.dll
