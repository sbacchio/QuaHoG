#ifndef _QHG_PROP_TRACE_H
#define _QHG_PROP_TRACE_H 1

static _Complex double
prop_trace(_Complex double X[NS*NC][NS*NC])
{
  _Complex double tr = 0;
  for(int cs=0; cs<NC*NS; cs++)
    tr += X[cs][cs];
      
  return tr;
}

#endif /* _QHG_PROP_TRACE_H */
