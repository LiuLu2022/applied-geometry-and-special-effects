#include "standardhidmanager.h"

// local
//#include "../glcontextsurfacewrapper.h"
#include "../gmlibwrapper.h"
#include "hidinputevent.h"
#include "hidkbmouseinput.h"

// gmlib
#include <gmSceneModule>
#include <gmParametricsModule>
using namespace GMlib;

// qt
#include <QtCore>
#include <QCoreApplication>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#define QT_NO_OPENGL
//#include <QtWidgets>

// stl
#include <cmath>
#include <cassert>

// Local Defines
#define SNAP 0.01



StandardHidManager::StandardHidManager( std::shared_ptr<GMlibWrapper> gmlib, QObject* parent )
  : HidManager(parent), _gmlib{gmlib} {

  // Mutable register
  _reg_wheel_state = false;


  // Default values
//  _previous_scene_pos = _current_scene_pos = Vector<int,2>(0.0f);

  _reg_next_key_event_type = KEY_NONE;
  _reg_next_mouse_event_type = MOUSE_NONE;
}

void StandardHidManager::registerKeyPressEvent(const QString& name, QKeyEvent *e) {

  registerRCPairName( name );
  registerKey( Qt::Key(e->key()), e->modifiers() );
  registerKeyEventType( KEY_PRESS );
  generateEvent();
}

void StandardHidManager::registerKeyReleaseEvent(const QString& name, QKeyEvent *e) {

  registerRCPairName( name );
  unregisterKey( Qt::Key(e->key()), e->modifiers() );
  registerKeyEventType( KEY_RELEASE );
  generateEvent();
}

void StandardHidManager::registerMouseDoubleClickEvent(const QString& name, QMouseEvent* e ) {

  registerRCPairName( name );
  registerWindowPosition( e->pos() );
  registerMouseButtons( e->buttons(), e->modifiers() );
  registerMouseEventType( MOUSE_DBL_CLICK );

  generateEvent();
}

void StandardHidManager::registerMouseMoveEvent(const QString& name, QMouseEvent* e) {

  registerRCPairName( name );
  registerWindowPosition( e->pos() );
  registerMouseEventType( MOUSE_MOVE );

  generateEvent();
}

void StandardHidManager::registerMousePressEvent(const QString& name, QMouseEvent* e) {

  registerRCPairName( name );
  registerWindowPosition( e->pos() );
  registerMouseButtons( e->buttons(), e->modifiers() );
  registerMouseEventType( MOUSE_CLICK );

  generateEvent();
}

void StandardHidManager::registerMouseReleaseEvent(const QString& name, QMouseEvent* e){

  registerRCPairName( name );
  registerWindowPosition( e->pos() );
  registerMouseButtons( e->buttons(), e->modifiers() );
  registerMouseEventType( MOUSE_RELEASE );

  generateEvent();
}


void StandardHidManager::registerKey(Qt::Key key, Qt::KeyboardModifiers modifiers) {

  _reg_keymods = modifiers;

  if( isKeyRegistered( key ) )
    return;

  _reg_keymap[key] = true;
}

void StandardHidManager::registerMouseButtons(Qt::MouseButtons buttons, Qt::KeyboardModifiers modifiers ) {

  _reg_mouse_buttons = buttons;
  _reg_keymods = modifiers;
}

void StandardHidManager::registerKeyEventType(StandardHidManager::KeyEventType type) {

  _reg_next_key_event_type = type;
}

void StandardHidManager::registerMouseEventType(StandardHidManager::MouseEventType type ) {

  _reg_next_mouse_event_type = type;
}

void StandardHidManager::setupDefaultHidBindings() {

  //// Register Hid Actions

  // Camera
//  QString ha_id_view_move_border_or_camera =
//      registerHidAction("View",
//                        "Move Border/Camera",
//                        "Move the border/camera. ",
//                        this, SLOT(heMoveBorderOrCamera()) );

  QString ha_id_view_move_camera =
      registerHidAction("View",
                        "Move Camera",
                        "Move the camera. "
                        "If not locked to the scene, it will pan the camera in the view plane. "
                        "If locked it will rotate the camera about the center of the scene." ,
                        this, SLOT(heMoveCamera(HidInputEvent::HidInputParams)) );

  QString ha_id_view_pan_h =
      registerHidAction("View",
                        "Pan Horizontally",
                        "Pan horizontally",
                        this, SLOT(hePanHorizontal(HidInputEvent::HidInputParams)) );

  QString ha_id_view_pan_v =
      registerHidAction("View",
                        "Pan Vertically",
                        "Pan vertically",
                        this, SLOT(hePanVertical(HidInputEvent::HidInputParams)) );

  QString ha_id_view_zoom=
      registerHidAction("View",
                        "Zoom",
                        "Zoom",
                        this, SLOT(heZoom(HidInputEvent::HidInputParams)) );

  QString ha_id_view_lock_to =
      registerHidAction("View",
                        "Lock To ...",
                        "Lock camera to an object or to the scene.",
                        this, SLOT(heLockTo(HidInputEvent::HidInputParams)) );

//  QString ha_id_view_moveborder =
//      registerHidAction( "View",
//                         "Move Border",
//                         "Move the border",
//                         this, SLOT(heMoveBorder()) );

//  QString ha_id_view_push_pop_viewset =
//      registerHidAction( "View",
//                         "Push/Pop View",
//                         "Push or pop the current view set.",
//                         this, SLOT(hePushPopViewSets()) );

  // Object Transformation
  QString ha_id_objtrans_scale =
      registerHidAction("Object transformation",
                        "Scale Objects",
                        "Scale objects",
                        this, SLOT(heScaleSelectedObjects(HidInputEvent::HidInputParams)) );

  QString ha_id_objtrans_move =
      registerHidAction("Object transformation",
                        "Move Objects",
                        "Move objects",
                        this, SLOT(heMoveSelectedObjects(HidInputEvent::HidInputParams)) );

  QString ha_id_objtrans_rotate =
      registerHidAction("Object transformation",
                        "Rotate Objects",
                        "Rotate objects",
                        this, SLOT(heRotateSelectedObjects(HidInputEvent::HidInputParams)) );

  // Object Selection
  QString ha_id_objsel_toggle_all =
      registerHidAction("Object selection",
                        "Toggle: Select all objects",
                        "Toggle Selection for all objects",
                        this, SLOT(heToggleSelectAllObjects()) );

  QString ha_id_objsel_select =
      registerHidAction("Object selection",
                        "Select Object",
                        "Select object",
                        this, SLOT(heSelectObject(HidInputEvent::HidInputParams)) );

  QString ha_id_objsel_select_multi =
      registerHidAction("Object selection",
                        "Select Objects",
                        "Select objects",
                        this, SLOT(heSelectObjects(HidInputEvent::HidInputParams)) );

  // Object Interaction
  QString ha_id_objint_toggle_edit =
      registerHidAction( "Object interaction",
                         "Toggle: Edit mode",
                         "Toggle edit mode for editable objects",
                         this, SLOT(heEdit()) );

  QString ha_id_objint_replot_high =
      registerHidAction( "Object interaction",
                         "Replot: QuickHigh",
                         "Replot with \"high\" resolution",
                         this, SLOT(heReplotQuickHigh()) );

  QString ha_id_objint_replot_med =
      registerHidAction( "Object interaction",
                         "Replot: QuickMedium",
                         "Replot with \"medium\" resolution",
                         this, SLOT(heReplotQuickMedium()) );

  QString ha_id_objint_replot_low =
      registerHidAction( "Object interaction",
                         "Replot: QuickLow",
                         "Replot with \"low\" resolution",
                         this, SLOT(heReplotQuickLow()) );



  // Rendering
  QString ha_id_render_toggle_shademode =
      registerHidAction( "Rendering",
                         "Toggle: Shading mode",
                         "Toggle shadeing mode",
                         this, SLOT(heToggleObjectDisplayMode()) );

  // Simulator
  QString ha_id_sim_toggle =
      registerHidAction( "Simulator",
                         "Toggle: Simulation",
                         "Toggle simulation",
                         this, SLOT(heToggleSimulation()) );

  // Various cleanup
  QString ha_id_var_lm_rel =
      registerHidAction( "Various",
                         "Left Mouse Release",
                         "Stuff that happens on left mouse release",
                         this, SLOT(heLeftMouseReleaseStuff()) );

  // Open/Close HidBindings view
  QString ha_id_var_open_close_hbview =
      registerHidAction( "Various",
                         "Open/Close Hid help",
                         "Toggle open/close the Hid bindings help view",
                         this, SLOT(heOpenCloseHidHelp()) );



  //// Set up initial mapping
  registerHidMapping( ha_id_objsel_toggle_all,            new KeyPressInput( Qt::Key_A ) );
  registerHidMapping( ha_id_objint_toggle_edit,           new KeyPressInput( Qt::Key_E ) );
  registerHidMapping( ha_id_objint_replot_high,           new KeyPressInput( Qt::Key_P, Qt::ShiftModifier ) );
  registerHidMapping( ha_id_objint_replot_med,            new KeyPressInput( Qt::Key_P ) );
  registerHidMapping( ha_id_objint_replot_low,            new KeyPressInput( Qt::Key_P, Qt::ControlModifier) );
  registerHidMapping( ha_id_sim_toggle,                   new KeyPressInput( Qt::Key_R ) );
  registerHidMapping( ha_id_render_toggle_shademode,      new KeyPressInput( Qt::Key_Z ) );

  registerHidMapping( ha_id_objsel_select,                new MousePressInput( Qt::RightButton ) );
  registerHidMapping( ha_id_view_lock_to,                 new MousePressInput( Qt::RightButton, Qt::ControlModifier ) );
  registerHidMapping( ha_id_objsel_select_multi,          new MousePressInput( Qt::RightButton, Qt::ShiftModifier ) );

//  registerHidMapping( ha_id_view_push_pop_viewset,        new MouseDoubleClickInput( Qt::LeftButton ) );

//  registerHidMapping( ha_id_view_move_border_or_camera,   new MouseMoveInput( Qt::LeftButton ) );
  registerHidMapping( ha_id_objtrans_scale,               new MouseMoveInput( Qt::LeftButton, Qt::ControlModifier | Qt::AltModifier ) );
  registerHidMapping( ha_id_objtrans_move,                new MouseMoveInput( Qt::LeftButton, Qt::ShiftModifier ) );
  registerHidMapping( ha_id_objtrans_rotate,              new MouseMoveInput( Qt::LeftButton, Qt::ControlModifier ) );

  registerHidMapping( ha_id_var_lm_rel,                   new MouseReleaseInput( Qt::LeftButton ) );
  registerHidMapping( ha_id_var_open_close_hbview,        new KeyPressInput( Qt::Key_Question, Qt::ShiftModifier ) );

  registerHidMapping( ha_id_view_pan_h,                   new WheelInput( Qt::ControlModifier ) );
  registerHidMapping( ha_id_view_pan_v,                   new WheelInput( Qt::ShiftModifier ) );
  registerHidMapping( ha_id_view_zoom,                    new WheelInput() );
}

GMlib::Point<int,2> StandardHidManager::toGMlibPoint(const QPoint& point) {

  return GMlib::Point<int,2>( point.x(), point.y() );
}

void StandardHidManager::unregisterKey(Qt::Key key, Qt::KeyboardModifiers modifiers) {


  // Update modifier
  _reg_keymods = modifiers;

  // Update keymap
  if( !isKeyRegistered(key) )
    return;

  _reg_keymap.remove( key );
}

void StandardHidManager::registerWheelData(bool state, int delta) {

  _reg_wheel_state = state;
  _reg_wheel_delta = delta;
}

void StandardHidManager::registerWheelEvent(const QString& name, QWheelEvent *e) {

  // Save position and wheel delta
  registerRCPairName( name );
  registerWindowPosition( e->pos() );
  registerWheelData( true, e->angleDelta().y() );

  generateEvent();
}

bool StandardHidManager::isKeyRegistered(Qt::Key key) const {

  return _reg_keymap.value(key, false);
}

bool StandardHidManager::isAnyKeysRegistered() const {

  return _reg_keymap.size() > 0;
}

bool StandardHidManager::isModKeyRegistered(Qt::KeyboardModifier keymod) const {

  // No modifier is a special case
  if( keymod == Qt::NoModifier ) {
    return _reg_keymods == Qt::NoModifier;
  }

  return (_reg_keymods & keymod) == keymod;
}

bool StandardHidManager::isMouseButtonRegistered(Qt::MouseButton button) const {

  // No button is a special case
  if( button == Qt::NoButton ) {
    return _reg_mouse_buttons == Qt::NoButton;
  }

  return (_reg_mouse_buttons & button) == button;
}

void StandardHidManager::generateEvent() {

  HidInputEvent::HidInputParams params;
  params["view_name"] = QVariant( _reg_rcpair_name );

  HidInputEvent::HidInputParams key_params {params};

  HidInputEvent::HidInputParams mouse_params {params};
  mouse_params["pos"]      = QVariant( _reg_view_pos );
  mouse_params["prev_pos"] = QVariant( _reg_view_prev_pos );

  HidInputEvent::HidInputParams wheel_params {params};
  wheel_params["wheel_delta"] = QVariant(_reg_wheel_delta);

  if( _reg_next_mouse_event_type == MOUSE_MOVE ) {
    QCoreApplication::sendEvent( this, new HidInputEvent( MouseMoveInput( _reg_mouse_buttons, _reg_keymods ), mouse_params ) );
    registerMouseEventType( MOUSE_NONE );
  }
  else if( _reg_wheel_state ) {
    QCoreApplication::sendEvent( this, new HidInputEvent( WheelInput( _reg_keymods ), wheel_params ) );
    registerWheelData(false,0);
  }
  else if( _reg_next_mouse_event_type != MOUSE_NONE ) {

    switch( _reg_next_mouse_event_type ) {
      case MOUSE_DBL_CLICK:
        QCoreApplication::sendEvent( this, new HidInputEvent( MouseDoubleClickInput( _reg_mouse_buttons, _reg_keymods ), mouse_params ) );
        break;
      case MOUSE_CLICK:
        QCoreApplication::sendEvent( this, new HidInputEvent( MousePressInput( _reg_mouse_buttons, _reg_keymods ), mouse_params ) );
        break;
      case MOUSE_RELEASE:
        QCoreApplication::sendEvent( this, new HidInputEvent( MouseReleaseInput( _reg_mouse_buttons, _reg_keymods ), mouse_params ) );
        break;
    }

    registerMouseEventType( MOUSE_NONE );
  }
  else if( _reg_next_key_event_type != KEY_NONE ) {

    switch( _reg_next_key_event_type ) {
      case KEY_PRESS: {
        QCoreApplication::sendEvent( this, new HidInputEvent( KeyPressInput( _reg_keymap, _reg_keymods ), key_params ) );
      } break;
      case KEY_RELEASE: {
        QCoreApplication::sendEvent( this, new HidInputEvent( KeyReleaseInput( _reg_keymap, _reg_keymods ), key_params ) );
      } break;
    }

    registerKeyEventType( KEY_NONE );
  }





//  qDebug() << "HidManager::processEvents()";
//  qDebug() << "  Keys:                  " << _no_keys_pressed;
//  qDebug() << "  Mouse buttons:         " << _mouse_buttons;
//  qDebug() << "  Modifiers:             " << _keymod;
//  qDebug() << "  Mouse move:            " << _mouse_move;
//  qDebug() << "  Mouse double click:    " << _mouse_dc;


  ///////////////////////
  // Mouse input handling
  ///////////////////////
//  if( getMouseMoveState() ) {

//    if( isMouseButtonPressed( Qt::LeftButton ) ) {
//      if( isModifierPressed( Qt::NoModifier ) ) {
//        if( _move_border )
//          heMoveBorder();
//        else
//          heMoveCamera();
//      }
//      else if( isModifierPressed( Qt::ControlModifier ) && isModifierPressed( Qt::AltModifier ) ) {
//        heScaleSelectedObjects();
//      }
//      else if( isModifierPressed( Qt::ShiftModifier ) ) {
//        heMoveSelectedObjects();
//      }
//      else if( isModifierPressed(Qt::ControlModifier ) ) {
//        heRotateSelectedObjects();
//      }
//    }
//    setMouseMoveState(false);
//  }
//  else if( getWheelState() ) {

//    // Pan "x" direction
//    if( isModifierPressed( Qt::ControlModifier ) )
//      hePanHorizontal();

//    // Pan "y" direction
//    else if( isModifierPressed( Qt::ShiftModifier ) )
//      hePanVertical();

//    // Zoom
//    else if( isModifierPressed( Qt::NoModifier ) )
//      heZoom();

//    setWheelState(false);
//  }
//  else {

//    if( getMouseDoubleClickState() && isMouseButtonPressed( Qt::LeftButton ) ) {

//      hePushPopViewSets();
//    }
//    else {

//      if( isMouseButtonPressed( Qt::LeftButton ) ) {

//        if( isModifierPressed( Qt::NoModifier ) ) {

//          Camera *cam = getActiveCamera();
//          if( !cam )
//            _move_border = true;
//        }
//      }
//      else if( isMouseButtonPressed( Qt::RightButton ) ) {

//        if( isModifierPressed( Qt::ShiftModifier ) ) {

//          heSelectObject();
//        }
//        else if( isModifierPressed( Qt::ControlModifier ) ) {

//          heLockCameraToObject();
//        }
//        else {

//          if( getSelectObject() )
//            heSelectObject();
//          else
//            heLockCameraToScene();
//        }
//      }
//      else {

//        if( _move_border )
//          _move_border = false;
//      }
//    }


//    /////////////////////
//    // Key input handling
//    /////////////////////
//    if( isKeysPressed() ) {

//      // Select/deselect all
//      if( isKeyPressed(Qt::Key_A ) )
//        heToggleSelectAllObjects();

//      // Edit
//      else if( isKeyPressed(Qt::Key_E) )
//        heEdit();

//      // Replot
//      else if( isKeyPressed(Qt::Key_P) ) {


//        if( isModifierPressed(Qt::ShiftModifier) )
//          heReplotQuickHigh();
//        else if( isModifierPressed(Qt::NoModifier) )
//          heReplotQuickLow();
//      }

//      // Quit
//      else if( isKeyPressed(Qt::Key_Q) )
//        heQuit();

//      // Toggle run
//      else if( isKeyPressed(Qt::Key_R) )
//        heToggleSimulation();

//      // Toggle shade mode for an object
//      else if( isKeyPressed(Qt::Key_Z) )
//        heToggleObjectDisplayMode();

//      // Grab framebuffer
//      else if( isKeyPressed(Qt::Key_F) )
//        emit signGrabFrameBuffer();
//    }
//  }

  // Save coords
//  savePos();
}






///////////////////////////////////////////////////////////////////////
//   *******   GMlib Scene Stuff
///////////////////////////////////////////////////////////////////////
void StandardHidManager::heDeSelectAllObjects() {

  scene()->removeSelections();
}

void StandardHidManager::heEdit() {

  const Array<SceneObject*> &sel_objs = scene()->getSelectedObjects();
  for( int i = 0; i < sel_objs.getSize(); i++ ) {

    SceneObject *sel_obj = sel_objs(i);

    // ERBS
    PERBSCurve<float> *ecObj = dynamic_cast<PERBSCurve<float>*>( sel_obj );
    PERBSSurf<float> *esObj = dynamic_cast<PERBSSurf<float>*>( sel_obj );
    PERBSTriangle<float> *etObj = dynamic_cast<PERBSTriangle<float>*>( sel_obj );

    // Bezier
    PBezierCurve<float> *bcObj = dynamic_cast<PBezierCurve<float>*>( sel_obj );
    PBezierSurf<float> *bsObj = dynamic_cast<PBezierSurf<float>*>( sel_obj );
    PBezierTriangle<float> *btObj = dynamic_cast<PBezierTriangle<float>*>( sel_obj );

    // Arc
    PArc<float> *acObj = dynamic_cast<PArc<float>*>( sel_obj );

    // ERBS
    if( ecObj ) {

      if( ecObj->isLocalCurvesVisible() )
        ecObj->hideLocalCurves();
      else
        ecObj->showLocalCurves();
    }
    else if( esObj ) {

      if( esObj->isLocalPatchesVisible() )
        esObj->hideLocalPatches();
      else
        esObj->showLocalPatches();
    }
    else if( etObj ) {
      if( etObj->isLocalPatchesVisible() )
        etObj->hideLocalPatches();
      else
        etObj->showLocalPatches();
    }
    // Bezier
    else if( bcObj ) {

      PERBSCurve<float> *parent = dynamic_cast<PERBSCurve<float>*>( bcObj->getParent() );
      if( parent ) {

        if( bcObj->toggleCollapsed() )
          bcObj->hideSelectors();
        else
          bcObj->showSelectors();
      }
      else {

        if( bcObj->isSelectorsVisible() )
          bcObj->hideSelectors();
        else
          bcObj->showSelectors();
      }
    }
    else if( bsObj ) {

      PERBSSurf<float> *parent = dynamic_cast<PERBSSurf<float>*>( bsObj->getParent());
      if( parent ) {

        if( bsObj->toggleCollapsed() )
          bsObj->hideSelectors();
        else
          bsObj->showSelectors(true);
      }
      else {

        if( bsObj->isSelectorsVisible() )
          bsObj->hideSelectors();
        else
          bsObj->showSelectors(true);
      }
    }
    else if( btObj ) {

      PERBSTriangle<float> *parent = dynamic_cast<PERBSTriangle<float>*>( btObj->getParent() );
      if( parent ) {

        if( btObj->toggleCollapsed() )
          btObj->hideSelectors();
        else
          btObj->showSelectors(true);
      }
      else {

        if( btObj->isSelectorsVisible() )
          btObj->hideSelectors();
        else
          btObj->showSelectors(true);
      }
    }
    else if( acObj ) {

      acObj->toggleCollapsed();
    }

  }
}

void StandardHidManager::heLockTo(const HidInputEvent::HidInputParams& params) {

  auto view_name = params["view_name"].toString();
  auto pos       = params["pos"].toPoint();

  auto cam     = findCamera(view_name);
  auto sel_obj = findSceneObject(view_name,pos);

  if( sel_obj )
    cam->lock( sel_obj );
  else if(cam->isLocked())
    cam->unLock();
  else {

    cam->lock(
      ( scene()->getSphereClean().getPos() - cam->getPos() ) *
      cam->getDir() );
  }
}

void StandardHidManager::heMoveCamera(const HidInputEvent::HidInputParams& params) {

  auto view_name = params["view_name"].toString();
  auto pos       = toGMlibPoint(params["pos"].toPoint());
  auto prev      = toGMlibPoint(params["prev_pos"].toPoint());

  auto *cam = findCamera(view_name);
  if( !cam )
    return;

  const float scale = cameraSpeedScale( cam );
  const Vector<float,2> delta (
     (pos(0) - prev(0)) * scale / cam->getViewportW(),
     (prev(1) - pos(1)) * scale / cam->getViewportH()
  );
  cam->move( delta );
}

void StandardHidManager::heMoveSelectedObjects( const HidInputEvent::HidInputParams& params ) {

  QString view_name = params["view_name"].toString();

  Camera *cam = findCamera(view_name);
  if( !cam )
    return;

  QPointF q_pos       = params["curr_pos"].toPointF();
  QPointF q_prev_pos  = params["prev_pos"].toPointF();

//  qDebug() << "//////////////////////";
//  qDebug() << "Pos: " << q_pos;
//  qDebug() << "Pos: " << q_prev_pos;

  const Vector<int,2> pos( q_pos.x(), q_pos.y() );
  const Vector<int,2> prev( q_prev_pos.x(), q_prev_pos.y() );
//  const Vector<int,2> pos = getPos();
//  const Vector<int,2> prev = getPPos();

  const Array<SceneObject*> &sel_objs = scene()->getSelectedObjects();
  for( int i = 0; i < sel_objs.getSize(); i++ ) {

    SceneObject* obj = sel_objs(i);

    const double dh = cam->deltaTranslate( obj );
    const Vector<float,3> deltav(
      ( ( prev(0) - pos(0) ) * dh ) * cam->getSide() +
      ( ( pos(1) - prev(1) ) * dh ) * cam->getUp() );

    if( deltav.getLength() > SNAP && deltav.getLength() < 1000.0 ) {

      if( obj->getTypeId() != GM_SO_TYPE_SELECTOR )
        obj->translateGlobal( deltav );
      else if( obj->getTypeId()== GM_SO_TYPE_SELECTOR )
        obj->editPos(deltav);
    }
  }
}

void StandardHidManager::hePanHorizontal(const HidInputEvent::HidInputParams& params) {

  QString view_name   = params["view_name"].toString();
  int     wheel_delta = params["wheel_delta"].toInt();

  Camera *cam = findCamera(view_name);
  if( cam )
    cam->move(
      Vector<float,2>(
        wheel_delta * cameraSpeedScale(cam) / cam->getViewportH(),
        0.0f
        ));
}

void StandardHidManager::hePanVertical(const HidInputEvent::HidInputParams& params) {

  QString view_name   = params["view_name"].toString();
  int     wheel_delta = params["wheel_delta"].toInt();

  Camera *cam = findCamera(view_name);
  if( cam )
    cam->move(
      Vector<float,2>(
        0.0f,
        wheel_delta * cameraSpeedScale(cam) / cam->getViewportW()
            ));
}

void StandardHidManager::heReplotQuick(int factor) {

  const Array<SceneObject*> &sel_objs = scene()->getSelectedObjects();

  for( int i = 0; i < sel_objs.getSize(); i++ ) {
    std::cout << "Selected object: " << sel_objs(i)->getIdentity() << std::endl;

    GMlib::SceneObject *sel_obj = sel_objs(i);

    GMlib::PCurve<float,3> *curve = dynamic_cast<GMlib::PCurve<float,3>*>( sel_obj );
    GMlib::PSurf<float,3> *surf = dynamic_cast<GMlib::PSurf<float,3>*>( sel_obj );

    if( curve ) {

      GMlib::PERBSCurve<float> *erbs = dynamic_cast<GMlib::PERBSCurve<float>*>(curve);
      if( erbs )
        erbs->replot(
          (erbs->getLocalCurves().getDim()-1)*factor + 1,
          1 );
      else
        curve->replot( std::pow<int>( factor, 2 ) * 100, 2 );
    }
    else if( surf ) {

      GMlib::PERBSSurf<float> *erbs = dynamic_cast<GMlib::PERBSSurf<float>*>(surf);
      if( erbs )
        erbs->replot(
          (erbs->getLocalPatches().getDim1()-1)*factor + 1,
          (erbs->getLocalPatches().getDim2()-1)*factor + 1,
          2, 2 );
      else {
        std::cout << "Replot psruf" << std::endl;
        surf->replot( 10 * factor, 10 * factor, 2, 2 );
      }
    }
  }
}

void StandardHidManager::heReplotQuickHigh() {

  heReplotQuick(20);
}

void StandardHidManager::heReplotQuickLow() {

  heReplotQuick(1);
}

void StandardHidManager::heReplotQuickMedium() {

  heReplotQuick(10);
}

void StandardHidManager::heRotateSelectedObjects(const HidInputEvent::HidInputParams& params) {

  auto view_name = params["view_name"].toString();
  auto pos       = toGMlibPoint(params["pos"].toPoint());
  auto prev      = toGMlibPoint(params["prev_pos"].toPoint());

  Camera *cam = findCamera(view_name);
  if( !cam )
    return;

  const Array<SceneObject*> &objs = scene()->getSelectedObjects();

  // Compute rotation axis and angle in respect to the camera and view.
  const UnitVector<float,3> rot_v =
    float( pos(0) - prev(0) ) * cam->getUp() -
    float( prev(1) - pos(1) ) * cam->getSide();
  const Angle ang(
    M_2PI * sqrt(
      pow( double( pos(0) - prev(0) ) / cam->getViewportW(), 2 ) +
      pow( double( prev(1) - pos(1) ) / cam->getViewportH(), 2 ) ) );


  int no_objs = 0;
  Sphere<float,3> sphere;
  for( int i = 0; i < objs.getSize(); ++i )
    if( objs(i)->getTypeId() != GM_SO_TYPE_SELECTOR ) {
      sphere += objs(i)->getSurroundingSphereClean();
      no_objs++;
    }

  for( int i = 0; i < objs.getSize(); ++i )
    if( objs(i)->getTypeId() != GM_SO_TYPE_SELECTOR )
      if( std::abs(pos(0)-prev(0)) > POS_TOLERANCE || std::abs(pos(1)-prev(1)) > POS_TOLERANCE )
        no_objs > 1 ? objs(i)->rotate( ang, sphere.getPos(), rot_v) : objs(i)->rotate( ang, rot_v);
}

void StandardHidManager::heScaleSelectedObjects(const HidInputEvent::HidInputParams& params) {

  auto view_name = params["view_name"].toString();
  auto pos       = toGMlibPoint(params["pos"].toPoint());
  auto prev      = toGMlibPoint(params["prev_pos"].toPoint());

  Camera *cam = findCamera(view_name);
  if( !cam )
    return;

  const Array<SceneObject*> &sel_objs = scene()->getSelectedObjects();
  for( int i = 0; i < sel_objs.getSize(); i++ ) {

    SceneObject* obj = sel_objs(i);

    const double dh = cam->deltaTranslate( obj );
    const Vector<float,3> deltav(
      ( ( prev(0) - pos(0) ) * dh ) * cam->getSide() +
      ( ( pos(1) - prev(1) ) * dh ) * cam->getUp() );


    if( deltav.getLength() < 1000.0 )
      obj->scale( Vector<float,3>( 1.0f + deltav(1) ) );
  }
}

void StandardHidManager::heSelectAllObjects() {

  Scene *scene = this->scene();
  for( int i = 0; i < scene->getSize(); ++i )
    heSelectObjectTree( (*scene)[i] );
}

void StandardHidManager::heSelectObject(const HidInputEvent::HidInputParams& params) {

  auto view_name = params["view_name"].toString();
  auto pos       = params["pos"].toPoint();

  auto obj = findSceneObject(view_name,pos);
  if( !obj )
    return;

  // Preserver object selection
  auto selected = obj->isSelected();
  heDeSelectAllObjects();
  obj->setSelected( !selected );
}

void StandardHidManager::heSelectObjects(const HidInputEvent::HidInputParams& params) {

  auto view_name = params["view_name"].toString();
  auto pos       = params["pos"].toPoint();

  if( auto obj = findSceneObject(view_name,pos) ) obj->toggleSelected();

//  if(obj) obj->toggleSelected();
}

void StandardHidManager::heSelectObjectTree( SceneObject* obj ) {

  // Do not select cameras or lights
  GMlib::Camera *cam   = dynamic_cast<GMlib::Camera*>( obj );
  GMlib::Light  *light = dynamic_cast<GMlib::Light*>( obj );
  if( !cam && !light )
    obj->setSelected(true);

  // Recursive Propagation
  for( int i = 0; i < obj->getChildren().getSize(); i++ )
    heSelectObjectTree( (obj->getChildren())[i] );
}

void StandardHidManager::heToggleObjectDisplayMode() {

  const Array<SceneObject*> &sel_objs = scene()->getSelectedObjects();

  qDebug() << "Toggling object display mode: " << sel_objs.getSize();
  for( int i = 0; i < sel_objs.getSize(); i++ ) {


    auto obj = sel_objs(i);
    GMlib::Array<GMlib::Visualizer*> &visus = obj->getVisualizers();
    for( int i = 0; i < visus.getSize(); i++ ) {

      qDebug() << "  obj: " << obj->getName() << " : " << reinterpret_cast<long int>(visus[i]);

      visus[i]->toggleDisplayMode();
    }
  }
}

void StandardHidManager::heToggleSimulation() {

  emit signToggleSimulation();
}

void StandardHidManager::heToggleSelectAllObjects() {

  if( scene()->getSelectedObjects().getSize() > 0 )
    heDeSelectAllObjects();
  else
    heSelectAllObjects();
}

//void StandardHidManager::heUnlockCamera() {

//  Camera *cam = activeCamera();
//  if( cam )
//    cam->unLock();
//}

void StandardHidManager::heZoom(const HidInputEvent::HidInputParams& params) {

  QString view_name       = params["view_name"].toString();
  int     wheel_delta     = params["wheel_delta"].toInt();

  // Qt comp scale
  wheel_delta /= 8;

  Camera *cam = findCamera(view_name);

  Camera *isocam = dynamic_cast<IsoCamera*>( cam );
  if( isocam ) {

    if( wheel_delta < 0 ) isocam->zoom( 1.05 );
    if( wheel_delta > 0 ) isocam->zoom( 0.95 );
  }
  else if( cam ) {

    double scale;
    if( cam->isLocked() )
      scale = cam->getLockDist();
    else
      scale = scene()->getSphere().getRadius();

    cam->move( wheel_delta*scale / cam->getViewportH() );
  }
}

void StandardHidManager::heLeftMouseReleaseStuff() {

  //  _move_border = false;
}

void StandardHidManager::heOpenCloseHidHelp() {

  emit signOpenCloseHidHelp();
  qDebug() << "Toggle Hid Help";
}

void StandardHidManager::registerRCPairName(const QString& name) {

  _reg_rcpair_name = name;
}

Camera* StandardHidManager::findCamera( const QString& view_name ) const {

  return _gmlib->getCamera(view_name.toStdString()).get();
}


float StandardHidManager::cameraSpeedScale(Camera *cam) const {

  if( !cam )
    return 1.0f;

  if(cam->isLocked())
    return M_2PI * cam->getLockDist();

  return scene()->getSphere().getRadius();
}



Scene* StandardHidManager::scene() const {

  return _gmlib->getScene().get();
}

SceneObject* StandardHidManager::findSceneObject( const QString& view_name, const QPoint& pos  ) {


  return _gmlib->findSceneObject( view_name, toGMlibViewPosition(view_name,pos) );
}

GMlib::Point<int,2> StandardHidManager::toGMlibViewPosition(const QString& view_name, const QPoint &pos) {

  auto cam = findCamera(view_name);
  return Vector<int,2>( int(pos.x()), cam->getViewportH() - int(pos.y()) - 1 );
}

void StandardHidManager::registerWindowPosition( const QPoint& pos ) {

  _reg_view_prev_pos = _reg_view_pos;
  _reg_view_pos = pos;
}
