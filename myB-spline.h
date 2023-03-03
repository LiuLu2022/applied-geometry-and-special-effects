#ifndef B_SMyBSpline_H
#define B_SMyBSpline_H

#include <parametrics/gmpcurve.h>


namespace GMlib {


  template <typename T>
  class MyBSpline : public PCurve<T,3> {
    GM_SCENEOBJECT(MyBSpline)

  public:


    MyBSpline( const DVector<Vector<T,3>>& c);        // SUB_CURVE type constructor
    MyBSpline( const DVector<Vector<T,3>>& p, int n ); // BEZIER_CURVE type constructor
    MyBSpline( const MyBSpline<T>& copy );
    virtual ~MyBSpline();



    //****************************************
    //****** Virtual public functions   ******
    //****************************************


    // from PCurve
    bool                   isClosed() const override;


  protected:


    DVector<Vector<T,3>>         _c;          //!< Local curves (control curves)
    DVector<T>                   _t;          //!< knot vector
    bool                         _cl;         //!< closed (or open) curve?
    int                          _d;
    int                          _k;




    // Virtual functions from PCurve, which have to be implemented locally
    void                   eval( T t, int d, bool l) const override;
    T                      getEndP()   const override;
    T                      getStartP() const override;

    // Local help functions
    Vector<T,3>           getB(int i, T t) const;
    T                     getW(int d, int i, T t) const;
    void                  knotvector(int n, int k);
    int                   findi(T t) const;


    //void                   init(bool closed);

  }; // END class MyBSMyBSpline



// Include MyBSMyBSpline class function implementations
//*****************************************
// Constructors and destructor           **
//*****************************************

/*! MyBSpline<T>::MyBSpline( const DVector<Vector<float,3>>& c  )
 *  Default constructor, to create a line starting in p and ending in p+v.
 *
 *  \param[in] p      The start point
 *  \param[in] v      A vector describing the line
 */
  template <typename T>
  inline
  MyBSpline<T>::MyBSpline( const DVector<Vector<T,3>>& c ) : PCurve<T,3>(20, 0, 0), _d(2), _k(3){
      _c = c;
//      _d = 2;
//      _k = _d+1;
      knotvector(c.getDim(), _k);
  }



  /*! MyBSpline<T>::MyBSpline( const Point<T,3>& p, const Vector<T,3>& v  )
   *  Default constructor, to create a line starting in p1 and ending in p2.
   *
   *  \param[in] p1      The start point
   *  \param[in] p2      The end point
   */
//  template <typename T>
//  inline
//  MyBSpline<T>::MyBSpline( const DVector<Vector<T,3>>& p, int n ) : PCurve<T,3>(20, 0, 0), _d(2), _k(3){

//      knotvector(n, _k)

//  }
  template <typename T>
  inline
  MyBSpline<T>::MyBSpline( const DVector<Vector<T,3>>& p, int n ) : PCurve<T,3>(20, 0, 0), _d(2), _k(3){

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
  //GMlib::DMatrix<T> A(c_points.getDim(),n);



  /*! MyBSpline<T>::MyBSpline(const MyBSpline<T>& copy )
   *  A copy constructor
   *  Making a copy of the curve (line)
   *
   *  \param[in] copy The curve to copy
   */
  template <typename T>
  inline
  MyBSpline<T>::MyBSpline( const MyBSpline<T>& copy ) : PCurve<T,3>(copy) {
      _d= copy._d;
      _k= copy._k;
      _t= copy._t;
      _c= copy._c;
  }



  /*! MyBSpline<T>::~MyBSpline()
   *  The destructor
   *  clean up and destroy all private data
   */
  template <typename T>
  MyBSpline<T>::~MyBSpline() {}




  //***************************************************
  // Overrided (public) virtual functons from PCurve **
  //***************************************************

  /*! bool MyBSpline<T>::isClosed() const
   *  To tell that this curve (line) is always open.
   *
   */
  template <typename T>
  bool MyBSpline<T>::isClosed() const {
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
  void MyBSpline<T>::eval( T t, int d, bool /*l*/ ) const {

    this->_p.setDim( d + 1 );

    int i = findi(t);
    Vector<T,3> b = getB(i, t);

    this -> _p[0] = b[0]*_c[i-2] + b[1]*_c[i-1] + b[2]*_c[i];

   // }
  }



  /*! T MyBSpline<T>::getStartP() const
   *  Provides the start parameter value associated with
   *  the eval() function implemented above.
   *  (the start parameter value = 0).
   *
   *  \return The parametervalue at start of the internal domain
   */
  template <typename T>
  T MyBSpline<T>::getStartP() const {
    return _t[_d];
  }



  /*! T MyBSpline<T>::getEndP() const
   *  Provides the end parameter value associated with
   *  the eval() function implemented above.
   *  (the end parameter value = 1).
   *
   *  \return The parametervalue at end of the internal domain
   */
  template <typename T>
  T MyBSpline<T>::getEndP() const {
    return _t[_c.getDim()];
  }

  template <typename T>
  void MyBSpline<T>::knotvector(int n, int k) {
      for(int i=0; i<_k; i++)
          _t.push_back(0);
      for(int i=_k; i<_c.getDim(); i++)
          _t.push_back(i-_d);
      for(int i=0; i<_k; i++)
          _t.push_back(_t[_c.getDim()-1]+1);

  }
  template <typename T>
  int MyBSpline<T>::findi(T t) const {
    for(int i=_d; i<_c.getDim(); i++)
          if(t>=_t[i] && t<_t[i+1])
          return i;
    return _c.getDim()-1;
} // END namepace GMlib

   template <typename T>
          T MyBSpline<T>::getW(int d, int i, T t) const{
          return (t-_t[i])/(_t[i+d]-_t[i]);
  }

   template <typename T>
   Vector<T,3> MyBSpline<T>::getB(int i, T t) const{

          T w1i = getW(1,i,t);
          T w2im1 = getW(2, i-1,t);
          T w2i = getW(2, i, t);
          Vector<T,3> b;
          b[0] = (1-w1i)*(1-w2im1);
          b[1] = (1-w1i)* w2im1 + w1i*(1-w2i);
          b[2] = w1i*w2i;
          return b;
  }

  }




#endif // B_SMyBSpline_H
