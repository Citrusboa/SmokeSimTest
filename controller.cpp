////////////////////////////////////////
// tester.cpp
////////////////////////////////////////

#include "core.h"

#define WINDOWTITLE	"Poly Drawer"

extern Controller* g_controller;

//set callback functions
static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}

static void resize(GLFWwindow *window, int x,int y)
{
	::g_controller->Resize(window, x, y);
}

static void keyboard(GLFWwindow * window, int key, int scancode,int action, int mods)		
{
	::g_controller->Keyboard(window, key, scancode, action, mods);
}

static void mousebutton(GLFWwindow *window,int button,int action,int mods)	
{
	::g_controller->MouseButton(window, button, action, mods);
}

static void mousemotion(GLFWwindow *window, double x, double y)
{
	::g_controller->MouseMotion(window, x, y);
}

void Controller::RegisterObject(Object* object)
{
	object->SetParentWindow(windowHandle_);
	object->Reset();
	objects_[objectNo_++] = object;
	activeObj_++;
}


Controller::Controller(int argc,char **argv) {
	winX_ = 640;
	winY_ = 480;

	activeObj_ = -1;
	objectNo_ = 0;
	maxObjectNo_ = 100;
	objects_ = new Object*[maxObjectNo_];

	isCtrlPressed_ = isLeftKeyPressed_ = isMiddleKeyPressed_ = isRightKeyPressed_ = false;
	isLightOn_ = true;
	prevX_ = prevY_ = 0;

	// Initialize components
	if (!glfwInit()) {
		cout << "glfwInit() failed!" << endl;
		exit(0);
	}

	// Create the window
	windowHandle_ = glfwCreateWindow(winX_, winY_, WINDOWTITLE, NULL, NULL);
	if (!windowHandle_) {
		cout << "Create Window failed"  << endl;
		exit(0);
	}
	glfwMakeContextCurrent(windowHandle_);
	glfwSetWindowPos(windowHandle_, 0, 0);

	// Background color
	glClearColor( 0., 0., 0., 1. );

	// Callbacks
	glfwSetErrorCallback(error_callback);
	glfwSetMouseButtonCallback(windowHandle_, mousebutton);
	glfwSetCursorPosCallback(windowHandle_, mousemotion);
	glfwSetKeyCallback(windowHandle_, keyboard);
	glfwSetWindowSizeCallback(windowHandle_, resize);

	InitCamera();
}

void Controller::InitCamera()
{
	camera_ = new Camera(windowHandle_);
}

void Controller::BeginLoop()
{
	while (!glfwWindowShouldClose(windowHandle_))
    {
        /* Render here */
		::g_controller->Render();

        /* Swap front and back buffers */
        glfwSwapBuffers(windowHandle_);

        /* Poll for and process events */
        glfwPollEvents();
    }
}

//////////////////////////////////////////////////////////////////////////////

Controller::~Controller() {
	glFinish();
	glfwDestroyWindow(windowHandle_);
	glfwTerminate();
}

void Controller::Reset() {
	camera_->Reset();

	for(int i = 0 ; i < objectNo_ ; i++)
		objects_[i]->Reset();
}

void Controller::Render() {
	//Begin drawing scene
	camera_->Reset();
	for(int i = 0 ; i < objectNo_ ; i++)
		objects_[i]->Show();

	//Finish drawing scene
	//glFinish();
}

int Controller::GetActiveObject(int mx, int my)
{
	//do some judgement;
	return activeObj_;
}

void Controller::Quit() {
	glFinish();
	glfwDestroyWindow(windowHandle_);
	exit(0);
}


void Controller::Resize(GLFWwindow *window, int x, int y) {
	winX_ = x;
	winY_ = y;
	
	for(int i = 0 ; i < objectNo_ ; i++)
		objects_[i]->Resize(window, x, y);
}

void Controller::Keyboard(GLFWwindow * window, int key, int scancode, int action, int mods)	
{
	if (action == GLFW_PRESS) {
		switch(key) {
			case GLFW_KEY_ESCAPE:		// Escape
				Quit();
				break;
			case GLFW_KEY_R:			//reset
				Reset();
				break;
			case GLFW_KEY_C:
				objects_[activeObj_]->ComputeValenceColor();
				/*
				isLightOn_ = !isLightOn_;
				if (isLightOn_) {
					glEnable(GL_LIGHT0);
					glEnable(GL_LIGHTING);
				}
				else {
					glDisable(GL_LIGHT0);
					glDisable(GL_LIGHTING);
				}
				*/
				break;
			case GLFW_KEY_LEFT_CONTROL:
				isCtrlPressed_ = true;
				break;
		}
	}
	else if(action == GLFW_RELEASE)
		switch (key) {
			case GLFW_KEY_LEFT_CONTROL:
				isCtrlPressed_ = false;
				break;
		}
}


void Controller::MouseButton(GLFWwindow *window, int button,int action,int mods) 
{

	//get active object and then transfer the message to the object
	if(action == GLFW_PRESS) {
		glfwGetCursorPos(window, &prevX_, &prevY_);
	}

	activeObj_ = GetActiveObject(prevX_, prevY_);
	objects_[activeObj_]->MouseButton(window, button, action, mods);
}


void Controller::MouseMotion(GLFWwindow *window, double nx, double ny) 
{
	objects_[activeObj_]->MouseMotion(window, nx, ny);
}

