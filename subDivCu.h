#ifndef CLS_SUBDIV_CURVE_H
#define CLS_SUBDIV_CURVE_H


#include <parametrics/gmpcurve.h>


namespace GMlib {


  template <typename T>
  class Cls_subdiv_curve : public PCurve<T,3> {
    GM_SCENEOBJECT(Cls_subdiv_curve)

  public:


    Cls_subdiv_curve( const DVector<Vector<T,3>>& c);        // SUB_CURVE type constructor
    //Cls_subdiv_curve( const DVector<Vector<T,3>>& p, int n ); // BEZIER_CURVE type constructor
    Cls_subdiv_curve( const Cls_subdiv_curve<T>& copy );
    virtual ~Cls_subdiv_curve();



    //****************************************
    //****** Virtual public functions   ******
    //****************************************


    // from PCurve
    bool                   isClosed() const override;
    void                   sample(int m, int d = 0);


  protected:


    DVector<Vector<T,3>>         _c;          //!< Local curves (control curves)
    //int                          _d;
    //int                          _k;




    // Virtual functions from PCurve, which have to be implemented locally
    void                   eval( T t, int d, bool l) const override;
    T                      getEndP()   const override;
    T                      getStartP() const override;

    // Local help functions
      void                  LaneRiesenfieldClosed(GMlib::DVector<GMlib::Vector<T,3>>& P,int k, int d);
      int                   double_part(std::vector<DVector<Vector<T,3>>>& P, int n);
      void                  smoothpartclosed(std::vector<DVector<Vector<T,3>>>& P,int k, int d);


//    Vector<T,3>           getB(int i, T t) const;
//    T                     getW(int d, int i, T t) const;
//    void                  knotvector(int n, int k);
//    int                   findi(T t) const;


    //void                   init(bool closed);

  }; // END class MyBSCls_subdiv_curve



// Include MyBSCls_subdiv_curve class function implementations
//*****************************************
// Constructors and destructor           **
//*****************************************

/*! Cls_subdiv_curve<T>::Cls_subdiv_curve( const DVector<Vector<float,3>>& c  )
 *  Default constructor, to create a line starting in p and ending in p+v.
 *
 *  \param[in] p      The start point
 *  \param[in] v      A vector describing the line
 */
  template <typename T>
  inline
  GMlib::Cls_subdiv_curve<T>::Cls_subdiv_curve( const DVector<Vector<T,3>>& c ) : PCurve<T,3>(20, 0, 0) {
      _c = c;
//      _d = 2;
//      _k = _d+1;
     // knotvector(c.getDim(), _k);
  }



  /*! Cls_subdiv_curve<T>::Cls_subdiv_curve( const Point<T,3>& p, const Vector<T,3>& v  )
   *  Default constructor, to create a line starting in p1 and ending in p2.
   *
   *  \param[in] p1      The start point
   *  \param[in] p2      The end point
   */
//  template <typename T>
//  inline
//  Cls_subdiv_curve<T>::Cls_subdiv_curve( const DVector<Vector<T,3>>& p, int n ) : PCurve<T,3>(20, 0, 0), _d(2), _k(3){

//      knotvector(n, _k)

//  }
  template <typename T>
  inline
  GMlib::Cls_subdiv_curve<T>::Cls_subdiv_curve( const Cls_subdiv_curve<T>& copy ) : PCurve<T,3>(20, 0, 0){

    _c= copy._c;
  }
  //GMlib::DMatrix<T> A(c_points.getDim(),n);



  /*! Cls_subdiv_curve<T>::Cls_subdiv_curve(const Cls_subdiv_curve<T>& copy )
   *  A copy constructor
   *  Making a copy of the curve (line)
   *
   *  \param[in] copy The curve to copy
   */
//  template <typename T>
//  inline
//  Cls_subdiv_curve<T>::Cls_subdiv_curve( const Cls_subdiv_curve<T>& copy ) : PCurve<T,3>(copy) {
//      _d= copy._d;
//      _k= copy._k;
//      _t= copy._t;
//      _c= copy._c;
//  }



  /*! Cls_subdiv_curve<T>::~Cls_subdiv_curve()
   *  The destructor
   *  clean up and destroy all private data
   */
  template <typename T>
  Cls_subdiv_curve<T>::~Cls_subdiv_curve() {}




  //***************************************************
  // Overrided (public) virtual functons from PCurve **
  //***************************************************

  /*! bool Cls_subdiv_curve<T>::isClosed() const
   *  To tell that this curve (line) is always open.
   *
   */
  template <typename T>
  bool Cls_subdiv_curve<T>::isClosed() const {
    return true;
  }



  //******************************************************
  // Overrided (protected) virtual functons from PCurve **
  //******************************************************

  /*! void PCircle<T>::eval( T t, int d, bool l ) const
   *  Evaluation of the curve at a given parameter value
   *  To compute position and d derivatives at parameter value t on the curve.
   *  7 derivatives are implemented
   *
   *  \param  t[in]  The parameter value to evaluate at
   *  \param  d[in]  The number of derivatives to compute
   *  \param  l[in]  (dummy) because left and right are always equal
   */
  template <typename T>
  void Cls_subdiv_curve<T>::eval( T t, int d, bool /*l*/ ) const {

    this->_p.setDim( d + 1 );

    }



  /*! T Cls_subdiv_curve<T>::getStartP() const
   *  Provides the start parameter value associated with
   *  the eval() function implemented above.
   *  (the start parameter value = 0).
   *
   *  \return The parametervalue at start of the internal domain
   */
  template <typename T>
  T Cls_subdiv_curve<T>::getStartP() const {
    return 0;
  }



  /*! T Cls_subdiv_curve<T>::getEndP() const
   *  Provides the end parameter value associated with
   *  the eval() function implemented above.
   *  (the end parameter value = 1).
   *
   *  \return The parametervalue at end of the internal domain
   */
  template <typename T>
  T Cls_subdiv_curve<T>::getEndP() const {
    return 1;
  }

  template <typename T>
   void GMlib::Cls_subdiv_curve<T>::LaneRiesenfieldClosed(GMlib::DVector<Vector<T,3>>& P,int k, int d) {

      int n = P.getDim();
      int m = pow(2,k)*n +1;
      std::cout << "n = " << n << std::endl;
      std::vector<DVector<Vector<T,3>>>& g= this->_visu[0].sample_val;
      g.resize(m);
      for(int i=0; i<n; i++){
          g[i][0]=P[i];
      }
          g[n++][0]=P[0][0];



      for(int i=0; i<k; i++){
          n=double_part(g,n);
          smoothpartclosed(g,n,d);
          std::cout << "n = " << n << std::endl;

      }
      Sphere<T, 3>& s=_visu[0].sur_sphere;
      s.reset();
      computeSurroundingSphere(g,s);




//      for(int i=0; i<m; i++)
//          s += g[i][0];
//      _visu[0].sur_sphere = s;
      //SceneObject::setSurroundingSphere(_visu[0].sur_sphere);

   }

//std::cout<< g << std::endl;

  }
  template <typename T>
  int GMlib::Cls_subdiv_curve<T>::double_part(std::vector<DVector<Vector<T,3>>>& P, int n) {

      for(int i=n-1; i>0; i--){
          P[2*i][0]=P[i][0];
          P[2*i-1][0] = (P[i][0]+P[i-1][0])/2;
      }
      return 2*n-1;

} // END namepace GMlib

   template <typename T>
   void GMlib::Cls_subdiv_curve<T>::smoothpartclosed(std::vector<DVector<Vector<T,3>>>& P,int n, int d){
          for(int j=1; j<d; j++){
              for(int i=0; i<n-1; i++){
                  P[i][0]= (P[i][0] + P[i+1][0])/2;

              }
              P[n-1][0]=P[0][0];

          }
  }

   template <typename T>
   void GMlib::Cls_subdiv_curve<T>::sample( int m, int d =0) {

     this-> _visu.resize(1);
     _checkSampleVal(m,d);
     //SceneObject::replot();
     LaneRiesenfieldClosed(_c, m, d);
     this->setEditDone();
   }






#endif // CLS_SUBDIV_CURVE_H
