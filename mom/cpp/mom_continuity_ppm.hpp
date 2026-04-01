// mom_continuity_ppm.hpp

#ifndef CONTINUITY_PPM_H_
#define CONTINUITY_PPM_H_

#include <AMReX_Box.H>

void ppm_limit_pos(Box const& bx,
                          Array4<Real> const& H_IN,
                          Array4<Real> & HL,
                          Array4<Real> & HR,
                          Real hmin);

void ppm_limit_cw84(Box const& bx,
                          Array4<Real> const& H_IN,
                          Array4<Real> & HL,
                          Array4<Real> & HR);

#endif
