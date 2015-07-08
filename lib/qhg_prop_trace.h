#ifndef _QHG_PROP_TRACE_H
#define _QHG_PROP_TRACE_H 1

static _Complex double
prop_trace(_Complex double X[NS*NC][NS*NC])
{
  _Complex double tr = 0.0;
  _Complex double temp;
  for(int cs=0; cs<NC*NS; cs++) {
    temp = tr + X[cs][cs];
    tr = temp;
  }
  
  return tr;
}

#endif /* _QHG_PROP_TRACE_H */
