#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include "nonlinfunc.hpp"
#include "timestepper.hpp"

using namespace ASC_ode;

#include "vector.hpp"
using namespace nanoblas;


template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  Vec<D> m_gravity=0.0;
public:
  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;
  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const override { return D*mss.masses().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override  
  {
    f = 0.0;

    auto xmat = x.asMatrix(mss.masses().size(), D);
    auto fmat = f.asMatrix(mss.masses().size(), D);

    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    for (auto spring : mss.springs())
      {
        auto [c1,c2] = spring.connectors;
        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        double force = spring.stiffness * (norm(p1-p2)-spring.length);
        Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);
        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += force*dir12;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= force*dir12;
      }

    for (size_t i = 0; i < mss.masses().size(); i++)
      fmat.row(i) *= 1.0/mss.masses()[i].mass;
  }

/*
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override   // Finite difference need to Modify !!!!!
  {
    // TODO: exact differentiation
    double eps = 1e-8;
    Vector<> xl(dimX()), xr(dimX()), fl(dimF()), fr(dimF());
    for (size_t i = 0; i < dimX(); i++)
      {
        xl = x;
        xl(i) -= eps;
        xr = x;
        xr(i) += eps;
        evaluate (xl, fl);
        evaluate (xr, fr);
        df.col(i) = 1/(2*eps) * (fr-fl);
      }
  }                               // Modify till here !!!!!
*/ 
    virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
        df = 0.0;

        auto xmat = x.asMatrix(mss.masses().size(), D);

        auto writeBlock = [&](size_t i, size_t j, const Mat<D,D>& M)
        {
            // 對 2D：row = 2*i + r, col = 2*j + c
            for (int r = 0; r < D; r++)
                for (int c = 0; c < D; c++)
                    df(i*D + r, j*D + c) += M(r,c);
        };

        for (auto &spring : mss.springs())
        {
            auto [c1, c2] = spring.connectors;

            // Position
            Vec<D> p1, p2;
            if (c1.type == Connector::FIX)
                p1 = mss.fixes()[c1.nr].pos;
            else
                p1 = xmat.row(c1.nr);

            if (c2.type == Connector::FIX)
                p2 = mss.fixes()[c2.nr].pos;
            else
                p2 = xmat.row(c2.nr);

            Vec<D> d = p2 - p1;
            double s = norm(d);
            if (s < 1e-12) continue;
            Vec<D> u = d / s;

            double k = spring.stiffness;
            double L = spring.length;

            Mat<D,D> A; // uu^T
            for (int r=0; r<D; r++)
                for (int c=0; c<D; c++)
                    A(r,c) = u[r]*u[c];

            Mat<D,D> I; // identity
            for (int r=0; r<D; r++)
                for (int c=0; c<D; c++)
                    I(r,c) = (r==c);

            Mat<D,D> B = I - A;

            // Derivative of a matrix（to p1 / p2）
            Mat<D,D> dF_dp1, dF_dp2;
            dF_dp1 = -k * ( A + ((s-L)/s) * B );
            dF_dp2 =  k * ( A + ((s-L)/s) * B );

            // divided by (accelerate = f = F/m)
            if (c1.type == Connector::MASS)
            {
                double m1 = mss.masses()[c1.nr].mass;
                Mat<D,D> M = (1.0/m1) * dF_dp1;
                writeBlock(c1.nr, c1.nr, M);

                if (c2.type == Connector::MASS)
                {
                    Mat<D,D> N = (1.0/m1) * (-dF_dp1);
                    writeBlock(c1.nr, c2.nr, N);
                }
            }

            if (c2.type == Connector::MASS)
            {
                double m2 = mss.masses()[c2.nr].mass;
                Mat<D,D> M = (1.0/m2) * dF_dp2;
                writeBlock(c2.nr, c2.nr, M);

                if (c1.type == Connector::MASS)
                {
                    Mat<D,D> N = (1.0/m2) * (-dF_dp2);
                    writeBlock(c2.nr, c1.nr, N);
                }
            }
        }
    }

};

#endif
