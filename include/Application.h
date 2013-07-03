#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Mesh.h"
#include "Real.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "DiscreteExteriorCalculus.h"
#include <iostream>

namespace DDG
{
   class Application
   {
   public:
       Application()
       {
           
       }
       
       void designVectorField(Mesh &mesh)
       {
           std::cout << std::endl;
           std::cout << "----------------------------------------" << std::endl;
           std::cout << "++++  start designing vector field  ++++" << std::endl;
           std::cout << "----------------------------------------" << std::endl;
           
           /* compute DEC operators */
           std::cout << std::endl;
           std::cout << "----computing DEC operators----" << std::endl;
           
           SparseMatrix<Real> star0, star1, star2;
           HodgeStar0Form<Real>::build(mesh, star0);
           HodgeStar1Form<Real>::build(mesh, star1);
           HodgeStar2Form<Real>::build(mesh, star2);
           
           SparseMatrix<Real> d0, d1;
           ExteriorDerivative0Form<Real>::build(mesh, d0);
           ExteriorDerivative1Form<Real>::build(mesh, d1);
           
           SparseMatrix<Real> del1 = star0.inverse() * d0.transpose() * star1;
           SparseMatrix<Real> del2 = star1.inverse() * d1.transpose() * star2;
           
           /* set up sink, source, and vortex constraints (no edge constraints yet) */
           DenseMatrix<Real> r_t = DenseMatrix<Real>(mesh.faces.size(), 1);
           r_t.zero(Real(0));
           DenseMatrix<Real> s_v = DenseMatrix<Real>(mesh.vertices.size(), 1);
           s_v.zero(Real(0));
           int numSources = 0, numSinks = 0, numVortices = 0;
           for (VertexIter v_it  = mesh.vertices.begin(); v_it != mesh.vertices.end(); v_it++)
           {
               if (v_it->sourceTag) {
                   numSources++;
                   s_v(v_it->index, 1) = -1000.f;
               }
               if (v_it->sinkTag) {
                   numSinks++;
                   s_v(v_it->index, 1) = 1000.f;
               }
           }
           for (FaceIter f_it = mesh.faces.begin(); f_it != mesh.faces.end(); f_it++)
           {
               if (f_it->vortexTag) {
                   numVortices++;
                   /* TODO: assign + or - based on orientation */
                   r_t(f_it->index, 1) = 1;
               }
           }
           for (VertexIter v_it  = mesh.vertices.begin(); v_it != mesh.vertices.end(); v_it++)
           {
               s_v(v_it->index, 1) = s_v(v_it->index, 1) - numSources*1000.f/(float)mesh.vertices.size() + numSinks*1000.f/(float)mesh.vertices.size();
           }
           for (FaceIter f_it = mesh.faces.begin(); f_it != mesh.faces.end(); f_it++)
           {
               if (f_it->vortexTag) {
                   r_t(f_it->index, 1) = r_t(f_it->index, 1) - numSinks*1.f/(float)mesh.faces.size();
               }
           }
           
           std::cout << "\tsources - " << numSources << std::endl;
           std::cout << "\tsinks - " << numSinks << std::endl;
           std::cout << "\tvortexes - " << numVortices << std::endl;
           
           /* compute M matrix */
           std::cout << std::endl;
           std::cout << "----computing giant M matrix----" << std::endl;
           
           SparseMatrix<Real> M = d1.transpose() * star2 * d1 + star1 * d0 * star0.inverse() * d0.transpose() * star1;
           
           /* compute edge coefficients */
           r_t = d1 * del2 * r_t;   // r_t: 2-form
           s_v = del1 * d0 * s_v;   // s_v: 0-form
           DenseMatrix<Real> ce, rhs;
           rhs = star1 * del2 * r_t + star1 * d0 * s_v;
           
           M = M + SparseMatrix<Real>::identity(mesh.edges.size());
           
           solve(M, ce, rhs);
           
           /* WHY IS THIS ALWAYS ZERO!?!?!? */
//           std::cout << "ce:" << std::endl;
//           for (int i = 0; i < ce.nRows(); i++) {
//               //std::cout << "\t| ";
//               for (int j = 0; j < ce.nColumns(); j++) {
//                   std::cout << ce(i, j) << ", ";
//               }
//               //std::cout << "|" << std::endl;
//           }
           
           /* assign vector field from edge coefficients */
           assignVectorField(mesh, ce);
           
           std::cout << std::endl;
           std::cout << "--------------------------------------" << std::endl;
           std::cout << "++++  end designing vector field  ++++" << std::endl;
           std::cout << "--------------------------------------" << std::endl;
       }
   private:
       void assignVectorField(Mesh &mesh, DenseMatrix<Real> ce)
       {
           /* parallel transport? */
           EdgeIter e = mesh.edges.begin();
           e->vec = e->he->flip->vertex->position - e->he->vertex->position;
           e++;
           for (e = mesh.edges.begin(); e != mesh.edges.end(); e++)
           {
               e->vec = e->he->flip->vertex->position - e->he->vertex->position;
               e->vec = ce(e->index, 1) * e->vec;
           }
       }
   };
}

#endif
