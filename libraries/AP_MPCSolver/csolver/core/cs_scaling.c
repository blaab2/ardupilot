#include "cs_scaling.h"

#include "../model/generated/turn_r5_probdata.h"

void cs_scaling_default(cs_scales *s)
{
    static const cs_real sx[6] = CS_PD_SCALE_X;
    static const cs_real su[2] = CS_PD_SCALE_U;
    int i;
    for (i = 0; i < 6; ++i)
        s->x[i] = sx[i];
    for (i = 0; i < 2; ++i)
        s->u[i] = su[i];
    s->T_ref = (cs_real)CS_PD_T_REF;
}
