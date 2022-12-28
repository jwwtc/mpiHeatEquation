#ifndef CONSTANTS_H
#define CONSTANTS_H

#define NX 129
#define NY 129

#define LX 20.0f
#define LY 20.0f

#define DX LX / ((REAL)(NX - 1.0))
#define DY LY / ((REAL)(NY - 1.0))

#define H DX
#define K   1.0f
#define DT  0.20f * DX * DX / K

#define TIME 200
#define STEPS (INT) ((double) (TIME) / (double) (DT))

#endif // CONSTANTS_H