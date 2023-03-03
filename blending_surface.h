#ifndef BLENDING_SURFACE_H
#define BLENDING_SURFACE_H

#include <parametrics/gmpcurve.h>
#include <parametrics/curves/gmpsubcurve.h>


#include "simple_sub_suf.h"



  template <typename T>
  class Blending_surface : public GMlib::PSurf<T,3> {
    GM_SCENEOBJECT(Blending_surface)
  public:
  // Constructors and destructor
    Blending_surface( GMlib::PSurf<T,3>* s, int nu, int nv);
    Blending_surface( const Blending_surface<T>& copy );

    virtual ~Blending_surface() {}

    bool                  isClosedU() const override;
    bool                  isClosedV() const override;

  protected:
    // Virtual functions from PSurf, which have to be implemented locally
    void          eval(T u, T v, int d1, int d2, bool lu = true, bool lv = true ) const override;
    T             getStartPU() const override;
    T             getEndPU()   const override;
    T             getStartPV() const override;
    T             getEndPV()   const override;




    // Protected data for the surface
    GMlib::PSurf<T,3>*                        _modelsurface;  // The original surface
    GMlib::DMatrix<GMlib::PSurf<T,3>*>        _c;          //!< Local curves (control curves)
    std::vector<T>                            _u;
    std::vector<T>                            _v; //!< knot vector
    int                                       _du;
    int                                       _ku;
    int                                       _dv;
    int                                       _kv;
    int                                       _nu;
    int                                       _nv;


  private:
    // Local help functions

    void                  knotvectoropen(std::vector<T>& _t, T st, T dt, T et, int n, int k); // making
    void                  knotvectorclosed(std::vector<T>& _t, T st, T dt, T et, int n, int k);
    GMlib::Vector<T,2>    getB(const std::vector<T>& _t, int i, T t) const; // others are using
    T                     getBfunc(T t) const;
    T                     getW(const std::vector<T>& _t,int d, int i, T t) const;
    GMlib::Vector<T,2>    getBd(const std::vector<T>& _t, int i, T t) const;
    T                     getBfuncd(T t) const;
    T                     getWd(const std::vector<T>& _t,int d, int i, T t) const;
    int                   findi(const std::vector<T>& _t, T t, bool cl) const;
    void                  localsurfaces(int nu, int nv, bool closed_u, bool closed_v);
    void                  localSimulate(double dt);

    // Private help functions
    void    set(GMlib::PSurf<T,3>* s, T su, T eu, T u, T sv, T ev, T v);

  }; // END class Blending_surface



  //***********************************************************
  // Include Blending_surface class function implementations ****
  //       sometimes located in a simplesubsurf.c file     ****
  //***********************************************************

  //*****************************************
  // Constructors and destructor           **
  //*****************************************

    template <typename T>
    inline
    Blending_surface<T>::Blending_surface( GMlib::PSurf<T,3>* s, int nu, int nv)
    {
        _modelsurface=s;  // The original surface
        _du=1;
        _ku=2;
        _dv=1;
        _kv=2;
        _nu=nu;
        _nv=nv;
        if(isClosedU())
            nu++;
        if(isClosedV())
            nv++;
        if (isClosedU())
            knotvectorclosed( _u, _modelsurface->getParStartU(), _modelsurface->getParDeltaU() , _modelsurface->getParEndU(), _nu, _ku);
        else
            knotvectoropen( _u, _modelsurface->getParStartU(), _modelsurface->getParDeltaU() , _modelsurface->getParEndU(), _nu, _ku);

        if (isClosedV())
            knotvectorclosed( _v, _modelsurface->getParStartV(), _modelsurface->getParDeltaV() , _modelsurface->getParEndV(), _nv, _kv);
        else
            knotvectoropen( _v, _modelsurface->getParStartV(), _modelsurface->getParDeltaV() , _modelsurface->getParEndV(), _nv, _kv);

        std::cout<< _u << " " << _v << std::endl;

        _c.setDim(nu,nv);

        for (int i=0; i<_nu; i++){
            for(int j=0; j<_nv; j++){
                _c[i][j]=new PSimpleSubSurf<T>(_modelsurface,_u[i],_u[i+2],_u[i+1],_v[j], _v[j+2], _v[j+1]);
                _c[i][j]->toggleDefaultVisualizer();
                _c[i][j]->sample(10,10,1,1);
                insert(_c[i][j]);
                _c[i][j]->setCollapsed(true);
            }

            if(isClosedV())
                _c[i][_nv] = _c[i][0];
        }
        if (isClosedU()){
            for(int j=0; j<_nv; j++)
                _c[_nu][j]= _c[0][j];
        }

        if (isClosedU() && isClosedV())
            _c[_nu][_nv]= _c[0][0];




    }





    template <typename T>
    inline
    Blending_surface<T>::Blending_surface( const Blending_surface<T>& copy ) : GMlib::PSurf<T,3>( copy )
    {
        _modelsurface = copy._modelsurface;  // The original surface
        _c = copy._c;          //!< Local curves (control curves)
        _u= copy._u;
        _v= copy._v; //!< knot vector
        _du= copy._du;
        _ku= copy._ku;
        _dv= copy._dv;
        _kv= copy._kv;
        _nu= copy._nu;
        _nv= copy._nv;

    }




    //*****************************************************
    // Overrided (protected) virtual functons from PSurf **
    //*****************************************************

    template <typename T>
    void Blending_surface<T>::eval( T u, T v, int d1, int d2, bool /*lu*/, bool /*lv*/) const {

      this->_p.setDim(1+d1, 1+d2);
      int i = findi(_u,u,isClosedU());
      int j = findi(_v,v,isClosedV());

      GMlib::Vector<T,2> bu = getB(_u, i, u);
      GMlib::Vector<T,2> bv = getB(_v, j, v);
      GMlib::Vector<T,2> bud = getBd(_u, i, u);//derivation of getb
      GMlib::Vector<T,2> bvd = getBd(_v, j, v);
      GMlib::DMatrix<GMlib::Vector<T,3>> c1 = _c(i-1)(j-1)->evaluateParent( u, v, d1, d2 );
      GMlib::DMatrix<GMlib::Vector<T,3>> c2 = _c(i)(j-1)->evaluateParent( u, v, d1, d2 );
      GMlib::DMatrix<GMlib::Vector<T,3>> c3 = _c(i-1)(j)->evaluateParent( u, v, d1, d2 );
      GMlib::DMatrix<GMlib::Vector<T,3>> c4 = _c(i)(j)->evaluateParent( u, v, d1, d2 );

      GMlib::Vector<T,3> M1= bu[0]*c1[0][0] + bu[1]*c2[0][0];
      GMlib::Vector<T,3> M2= bu[0]*c3[0][0] + bu[1]*c4[0][0];
      GMlib::Vector<T,3> M1u= bud[0]*c1[0][0] + bu[0]*c1[0][1] + bud[1]*c2[0][0] + bu[1]*c2[0][1]; // c1[0][1] = c1u=c2u=  [0][1]
      GMlib::Vector<T,3> M2u= bud[0]*c3[0][0] + bu[0]*c3[0][1] + bud[1]*c4[0][0] + bu[1]*c4[0][1];
      GMlib::Vector<T,3> M1v= bu[0]*c1[1][0] + bu[1]*c2[1][0];
      GMlib::Vector<T,3> M2v= bu[0]*c3[1][0] + bu[1]*c4[1][0];

      this->_p[0][0] = bv[0]*M1 +bv[1]*M2;
      this->_p[0][1] = bv[0]*M1u + bv[1]*M2u;
      this->_p[1][0] = bvd[0]*M1 + bv[0]*M1v + bvd[1]*M2 +bv[1]*M2v;

    }


    template <typename T>
    T Blending_surface<T>::getStartPU() const {
      return _modelsurface->getParStartU();
    }


    template <typename T>
    T Blending_surface<T>::getEndPU() const {
      return _modelsurface->getParEndU();
    }


    template <typename T>
    T Blending_surface<T>::getStartPV() const {
      return _modelsurface->getParStartV();
    }


    template <typename T>
    T Blending_surface<T>::getEndPV() const {
      return _modelsurface->getParEndV();
    }



    //***************************
    // Private help functions  **
    //***************************

    template <typename T>
    inline
    void Blending_surface<T>::set(GMlib::PSurf<T,3>* s, T su, T eu, T u, T sv, T ev, T v) {
      _s  = s;
      _su = su;
      _sv = sv;
      _eu = eu;
      _ev = ev;
      _u  = u;
      _v  = v;
    }

    template <typename T>
    bool Blending_surface<T>::isClosedU() const {
      return _modelsurface->isClosedU();
    }
    template <typename T>
    bool Blending_surface<T>::isClosedV() const {
      return _modelsurface->isClosedV();
    }


    template <typename T>
    GMlib::Vector<T,2> Blending_surface<T>::getB(const std::vector<T>& _t, int i, T t) const{

           T w1i = getBfunc(getW(_t, 1, i, t));
           GMlib::Vector<T,2> b;
           b[0] = (1-w1i);
           b[1] = w1i;
           return b;
   }

    template <typename T>
    T Blending_surface<T>::getBfunc(T t) const{
        //return 3*t*t-2*t*t*t;
        return t-(1/(2*M_PI))*sin(2*M_PI*t);//2nd order trigonometric function
    }

    template <typename T>
    T Blending_surface<T>::getW(const std::vector<T>& _t, int d, int i, T t) const{
        T w=(t-_t[i])/(_t[i+d]-_t[i]);
        return w;
    }

    template <typename T>
    T Blending_surface<T>::getWd(const std::vector<T>& _t, int d, int i, T t) const{
        T wd=1/(_t[i+d]-_t[i]);
        return getBfuncd(getW(_t,d,i,t))*wd;
    }

    template <typename T>
    T Blending_surface<T>::getBfuncd(T t) const{
        return 1-cos(2*M_PI*t);//2nd order trigonometric function
    }

    template <typename T>
    GMlib::Vector<T,2> Blending_surface<T>::getBd(const std::vector<T>& _t, int i, T t) const{

           T w1i = getWd(_t, 1, i, t);
           GMlib::Vector<T,2> b;
           b[0] = -w1i;
           b[1] = w1i;
           return b;
   }

    template <typename T>
    int Blending_surface<T>::findi(const std::vector<T>& _t, T t, bool cl) const {

        int n=_t.size()-2;
        if (cl)
            n++;
        for(int i=_du; i<n; i++)
            if(t>=_t[i] && t<_t[i+1])
                return i;
        return n-1;
    }

    template <typename T>
    void Blending_surface<T>::knotvectoropen(std::vector<T>& _t, T st, T dt, T et, int n, int k) {

        dt /= (n-1); // small pieces

        for(int i=0; i<k; i++)
            _t.push_back(st);
        for(int i=k; i<n; i++)
            _t.push_back(st+(i-1)*dt);
        for(int i=0; i<k; i++)
            _t.push_back(et);

    }

    template <typename T>
    void Blending_surface<T>::knotvectorclosed(std::vector<T>& _t, T st, T dt, T et, int n, int k) {

        dt /= n;

        for(int i=0; i<k; i++)
            _t.push_back(st);
        for(int i=k; i<n; i++)
            _t.push_back(st+(i-1)*dt);
        for(int i=0; i<k; i++)
            _t.push_back(et);
        _t[0]=_t[1]-dt;
        _t[n+1]=_t[n]+dt;
    }

    template <typename T>
    void Blending_surface<T>::localsurfaces(int nu, int nv, bool closed_u, bool closed_v){

        _c.setDim(nu + isClosedU()?1:0, nv + isClosedV()?1:0);
        for( int j=0; j<nu; j++){
            for(int i=0; i<nv; i++){
                PSimpleSubSurf<T>* ps = new PSimpleSubSurf<T>(_modelsurface, _u[i],_u[i+2],_u[i+1], _v[j],_v[j+2],_v[j+1]);
                _c[i][j]=ps;
                ps->toggleDefaultVisualizer();
                ps->sample(10,10, 1, 1);
                this->insert(ps);
                ps->setCollapsed(true);
            }

            if(isClosedU())
                _c[nu][j]=_c[0][j];

        }
         if (isClosedV()){
             for(int i=0; i<nv; i++)
                 _c[i][nv]=_c[i][nv];
         }
         else if(isClosedU() && isClosedV())
             _c[nu][nv]=_c[0][0];

    }



    template <typename T>
    void Blending_surface<T>::localSimulate(double dt){




          sample(20,20,1,1);
    }







//#endif // GM_PARAMETRICS_SURFACES_Blending_surface_H

#endif // BLENDING_SURFACE_H
