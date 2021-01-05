#ifndef OpenGL
#include<Gl/glew.h>
#include<GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
using namespace std;

class OpenGL {
    private:
        GLFWwindow* window;
        const int xpos = 299;
        const int ypos = 20;
        GLuint vertexshader;
        GLuint fragmentshader;
        GLchar infolog[300];
        GLuint shaderProgram;
        GLint s_Compile;
        GLuint Node_buffer;
        GLuint VAO;
        static const GLchar *vertex_shader;
        static const GLchar *fragment_shader;
    public:
        OpenGL();                                    /* Initialize OpenGL libraries */
        void poll();                                 /* Poll clear screen */
        void post();                                 /* Swap buffers */
        int close();                                 /* Close window */
        void vshader();                              /* Compile vertex shader */
        void fshader();                              /* Compile fragment shader */
        void linkShaders();                          /* Link shaders */
        void genBuffer();                            /* Generate Buffers */
        void delBuffer();                            /* Delete Buffers*/
        void render(std::vector<GLfloat> Vertices);
};

// 2d vector
struct vec2{
    float *a;
    float *b;
    float x(){return *a;}
    float y(){return *b;}
};

// Some function declarations
void key_callback(GLFWwindow*window, int key, int scancode, int action, int mode);
float dot(vec2,vec2);

#endif
