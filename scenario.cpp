
#include <iostream>

#include "scenario.h"
#include "testtorus.h"

#include <parametrics/curves/gmpline.h>
#include <parametrics/surfaces/gmpplane.h>
#include <parametrics/curves/gmpcircle.h>
#include "parametric_curve.h"
#include "myB-spline.h"
#include "circle.h"
#include "subDivCu.h"
#include "blending_spline_curve.h"
#include "blending_surface.h"
#include "simple_sub_suf.h"

// hidmanager
#include "hidmanager/defaulthidmanager.h"

// gmlib
#include <scene/light/gmpointlight.h>
#include <scene/sceneobjects/gmpathtrack.h>
#include <scene/sceneobjects/gmpathtrackarrows.h>

// qt
#include <QQuickItem>


template <typename T>
inline
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  out << v.size() << std::endl;
  for(uint i=0; i<v.size(); i++) out << " " << v[i];
  out << std::endl;
  return out;
}




void Scenario::initializeScenario() {

  // Insert a light
  GMlib::Point<GLfloat,3> init_light_pos( 2.0, 4.0, 10 );
  GMlib::PointLight *light = new GMlib::PointLight(  GMlib::GMcolor::white(), GMlib::GMcolor::white(),
                                                     GMlib::GMcolor::white(), init_light_pos );
  light->setAttenuation(0.8f, 0.002f, 0.0008f);
  this->scene()->insertLight( light, false );

  // Insert Sun
  this->scene()->insertSun();

  // Default camera parameters
  int init_viewport_size = 600;
  GMlib::Point<float,3>  init_cam_pos( 0.0f, 0.0f, 0.0f );
  GMlib::Vector<float,3> init_cam_dir( 0.0f, 1.0f, 0.0f );
  GMlib::Vector<float,3> init_cam_up(  1.0f, 0.0f, 0.0f );

  // Projection cam
  auto proj_rcpair = createRCPair("Projection");
  proj_rcpair.camera->set(init_cam_pos,init_cam_dir,init_cam_up);
  proj_rcpair.camera->setCuttingPlanes( 1.0f, 8000.0f );
  proj_rcpair.camera->rotateGlobal( GMlib::Angle(-45), GMlib::Vector<float,3>( 1.0f, 0.0f, 0.0f ) );
  proj_rcpair.camera->translateGlobal( GMlib::Vector<float,3>( 0.0f, -20.0f, 20.0f ) );
  scene()->insertCamera( proj_rcpair.camera.get() );
  proj_rcpair.renderer->reshape( GMlib::Vector<int,2>(init_viewport_size, init_viewport_size) );


  /***************************************************************************
   *                                                                         *
   * Standar example, including path track and path track arrows             *
   *                                                                         *
   ***************************************************************************/

  GMlib::Material mm(GMlib::GMmaterial::polishedBronze());
  mm.set(45.0);

//  auto ptom = new TestTorus(1.0f, 0.4f, 0.6f);
//  ptom->toggleDefaultVisualizer();
//  ptom->sample(60,60,1,1);
//  this->scene()->insert(ptom);
//  auto ptrack = new GMlib::PathTrack();
//  ptrack->setLineWidth(2);
//  ptom->insert(ptrack);
//  auto ptrack2 = new GMlib::PathTrackArrows();
//  ptrack2->setArrowLength(2);
//  ptom->insert(ptrack2);

 // -----------parametric curve--------
  auto ptom = new GMlib::parametric_curve<float>();
  ptom->toggleDefaultVisualizer();
  ptom->sample(160,1);
  this->scene()->insert(ptom);

  //------------B-spline------------
//  GMlib::DVector<GMlib::Vector<float, 3>> P(50);
//  auto cir = new GMlib::PCircle<float>(10);
//  for (int i = 0; i < 50; ++i) {
//      P[i]= cir->getPosition(cir->getParStart()+i* cir->getParDelta()/49);
//  }
//  auto mybspl = new GMlib::MyBSpline<float>(P,7);
//  mybspl->toggleDefaultVisualizer();
//  mybspl->sample(50,0);
//  mybspl->setColor(GMlib::GMcolor::azure());
//  mybspl->setLineWidth(5);
//  this->scene()->insert(mybspl);

  //------------circle------------
//  GMlib::DVector<GMlib::Vector<float, 3>> P(50);
//  auto cir = new GMlib::PCircle<float>(10);
//  for (int i = 0; i < 50; ++i) {
//      P[i]= cir->getPosition(cir->getParStart()+i* cir->getParDelta()/49);
//  }
//  auto mybspl = new GMlib::circle<float>(P,7);
//  mybspl->toggleDefaultVisualizer();
//  mybspl->sample(60,0);
//  mybspl->setColor(GMlib::GMcolor::azure());
//  mybspl->setLineWidth(5);
//  this->scene()->insert(mybspl);

  //----------Sub division curve-----------
//  GMlib::DVector<GMlib::Vector<float, 3>> c_points(6);
//   c_points[0] = GMlib::Vector<float,3>(0,0,0);
//   c_points[1] = GMlib::Vector<float,3>(0,2,0);
//   c_points[2] = GMlib::Vector<float,3>(1,2,0);
//   c_points[3] = GMlib::Vector<float,3>(4,2,0);
//   c_points[4] = GMlib::Vector<float,3>(2,0,0);
//   c_points[5] = GMlib::Vector<float,3>(5,1,0);
//    auto ptom = new GMlib::Cls_subdiv_curve<float>(c_points);
//    ptom->toggleDefaultVisualizer();
//    ptom->sample(3,2);
//    ptom->setColor(GMlib::GMcolor::yellow());
//    ptom->setLineWidth(5);
//    this->scene()->insert(ptom);

  //----------blending spline curve-----------
//    auto ptom = new GMlib::parametric_curve<float>();
//      ptom->toggleDefaultVisualizer();
//      ptom->sample(80,1);
//      ptom->setColor(GMlib::GMcolor::azure());
//      ptom->setLineWidth(5);
//      this->scene()->insert(ptom);

//      auto bscurve = new GMlib::Blending_curve<float>(ptom,4);
//      bscurve->toggleDefaultVisualizer();
//      bscurve->sample(80,1);
//      bscurve->setColor(GMlib::GMcolor::azure());
//      bscurve->setLineWidth(5);
//      this->scene()->insert(bscurve);

  //----------blending surface-----------
//      GMlib::Point<float, 3> myPoint(0, 0, 0);
//        GMlib::Vector<float, 3> V1(1, 0, 0);
//        GMlib::Vector<float, 3> V2(0, 1, 0);
//        auto myPlane = new GMlib::PPlane<float>(myPoint, V1, V2);
//        myPlane->toggleDefaultVisualizer();
//        myPlane->sample(20,20,1,1);
//        //this->scene()->insert(myPlane);
//        auto myBlendingSurface = new Blending_surface<float>(myPlane, 3,3);
//        myBlendingSurface->toggleDefaultVisualizer();
//        myBlendingSurface->sample(20,20,1,1);
//        this->scene()->insert(myBlendingSurface);

}




void Scenario::cleanupScenario() {

}




void Scenario::callDefferedGL() {

  GMlib::Array< const GMlib::SceneObject*> e_obj;
  this->scene()->getEditedObjects(e_obj);

  for(int i=0; i < e_obj.getSize(); i++)
    if(e_obj(i)->isVisible()) e_obj[i]->replot();
}

