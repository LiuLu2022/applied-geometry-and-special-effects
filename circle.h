#ifndef CIRCLE_H
#define CIRCLE_H
#include <parametrics/gmpcurve.h>


namespace GMlib {


  template <typename T>
  class circle : public PCurve<T,3> {
    GM_SCENEOBJECT(circle)
  public:
    circle();
    circle( const circle<T>& copy );
    virtual ~circle();

    //****************************************
    //****** Virtual public functions   ******
    //****************************************

    // from PCurve
    bool                isClosed() const override;

  protected:
    // Virtual functions from PCurve, which have to be implemented locally
    void                eval(T t, int d, bool l) const override;
    T                   getStartP() const override;
    T                   getEndP()   const override;

    // Protected data for the curve
    Point<T,3>      _pt;
    Vector<T,3>     _v;

  }; // END class circle



// Include circle class function implementations

//*****************************************
// Constructors and destructor           **
//*****************************************


  /*! circle<T>::circle( const Point<T,3>& p, const Vector<T,3>& v  )
   *  Default constructor, to create a line starting in p1 and ending in p2.
   *
   *  \param[in] p1      The start point
   *  \param[in] p2      The end point
   */
  template <typename T>
  inline
  GMlib::circle<T>::circle() : PCurve<T,3>(20, 0, 1){}



  /*! circle<T>::circle(const circle<T>& copy )
   *  A copy constructor
   *  Making a copy of the curve (line)
   *
   *  \param[in] copy The curve to copy
   */
  template <typename T>
  inline
  circle<T>::circle( const circle<T>& copy ) : PCurve<T,3>(copy) {
      _c.setDim(n);
           knotvector(n, _k);
           //std::cout<<_t<<std::endl;
           int m = p.getDim();
           DMatrix<T>A(m,n,T(0));
           for(int j=0; j<m;j++){
               T t =getParStart()+j*getParDelta()/(m-1);
               int i = findi(t);
               //std::cout<<t<< " "<<i<<std::endl;
               GMlib::Vector<T,3> b=getB(i,t);
               A[j][i-2]=b[0];
               A[j][i-1]=b[1];
               A[j][i]=b[2];
           }
           //std::cout<<A<<std::endl;
           GMlib::DMatrix<T>             AT=A;
           AT.transpose();
           GMlib::DMatrix<T>             B=AT*A;
           GMlib::DVector<Vector<T,3>>   q=AT*p;
           B.invert();
           _c=B*q;



  }



  /*! circle<T>::~circle()
   *  The destructor
   *  clean up and destroy all private data
   */
  template <typename T>
  circle<T>::~circle() {}




  //***************************************************
  // Overrided (public) virtual functons from PCurve **
  //***************************************************

  /*! bool circle<T>::isClosed() const
   *  To tell that this curve (line) is always open.
   *
   */
  template <typename T>
  bool circle<T>::isClosed() const {
    return false;
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
  void circle<T>::eval( T t, int d, bool /*l*/ ) const {

    this->_p.setDim( d + 1 );

    this->_p[0][0] = 6*t -9*t*t +6*t*t*t ;
    this->_p[0][1] = -6*t + 18*t*t -12*t*t*t;
    this->_p[0][2] = 0;

    if( this->_dm == GM_DERIVATION_EXPLICIT ) {
      if( d )     this->_p[1] = _v;
      if( d > 1 ) this->_p[2] = Vector<T,3>(T(0));
      if( d > 2 ) this->_p[3] = Vector<T,3>(T(0));
      if( d > 3 ) this->_p[4] = Vector<T,3>(T(0));
      if( d > 4 ) this->_p[5] = Vector<T,3>(T(0));
      if( d > 5 ) this->_p[6] = Vector<T,3>(T(0));
      if( d > 6 ) this->_p[7] = Vector<T,3>(T(0));
    }
  }



  /*! T circle<T>::getStartP() const
   *  Provides the start parameter value associated with
   *  the eval() function implemented above.
   *  (the start parameter value = 0).
   *
   *  \return The parametervalue at start of the internal domain
   */
  template <typename T>
  T circle<T>::getStartP() const {
    return T(0);
  }



  /*! T circle<T>::getEndP() const
   *  Provides the end parameter value associated with
   *  the eval() function implemented above.
   *  (the end parameter value = 1).
   *
   *  \return The parametervalue at end of the internal domain
   */
  template <typename T>
  T circle<T>::getEndP() const {
    return T(1);
  }
} // END namepace GMlib
#endif // CIRCLE_H
