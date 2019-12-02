{
   const double tmp0 = (1.0 / (ADM_3METRIC[LAPSE])*(ADM_3METRIC[LAPSE]));
   const double tmp1 = ADM_3METRIC[SHIFTX]*tmp0;
   const double tmp2 = ADM_3METRIC[SHIFTY]*tmp0;
   const double tmp3 = ADM_3METRIC[SHIFTZ]*tmp0;
   const double tmp4 = ADM_3METRIC[GUPXY] - ADM_3METRIC[SHIFTY]*tmp1;
   const double tmp5 = ADM_3METRIC[GUPXZ] - ADM_3METRIC[SHIFTZ]*tmp1;
   const double tmp6 = ADM_3METRIC[GUPYZ] - ADM_3METRIC[SHIFTZ]*tmp2;
   g4up[0][0] = -tmp0;
   g4up[0][1] = tmp1;
   g4up[0][2] = tmp2;
   g4up[0][3] = tmp3;
   g4up[1][0] = tmp1;
   g4up[1][1] = ADM_3METRIC[GUPXX] - (ADM_3METRIC[SHIFTX])*(ADM_3METRIC[SHIFTX])*tmp0;
   g4up[1][2] = tmp4;
   g4up[1][3] = tmp5;
   g4up[2][0] = tmp2;
   g4up[2][1] = tmp4;
   g4up[2][2] = ADM_3METRIC[GUPYY] - (ADM_3METRIC[SHIFTY])*(ADM_3METRIC[SHIFTY])*tmp0;
   g4up[2][3] = tmp6;
   g4up[3][0] = tmp3;
   g4up[3][1] = tmp5;
   g4up[3][2] = tmp6;
   g4up[3][3] = ADM_3METRIC[GUPZZ] - (ADM_3METRIC[SHIFTZ])*(ADM_3METRIC[SHIFTZ])*tmp0;
}
