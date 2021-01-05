#include "OpenGL.h"
OpenGL::OpenGL(){
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
    glfwWindowHint(GLFW_OPENGL_PROFILE,GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE,GLFW_FALSE);
    glfwWindowHint(GLFW_SAMPLES,4);
    glEnable(GL_MULTISAMPLE);

    window = glfwCreateWindow(768,768,"A Most Excellent Window",0,0);
    if(!window){glfwTerminate();cout<<"Could not create window\n";}
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if(glewInit()!=GLEW_OK){cout<<"Could not initialize GLEW\n";}

    glfwSetInputMode(window,GLFW_STICKY_KEYS,GL_TRUE);
    glViewport(0,0,768,768);
    glfwSetWindowPos(window,xpos,ypos);

    glfwSetKeyCallback(window,key_callback);

    vshader();
    fshader();
    linkShaders();
}
void OpenGL::poll(){
    glfwPollEvents();
    glClearColor(0.5f, 0.5f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}
void OpenGL::post(){
    glfwSwapBuffers(window);
}
int OpenGL::close(){
    return glfwWindowShouldClose(window);
}
void OpenGL::vshader(){
    vertexshader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexshader,1,&vertex_shader,0);
    glCompileShader(vertexshader);;
    glGetShaderiv(vertexshader, GL_COMPILE_STATUS,&s_Compile);
    if(!s_Compile){
        glGetShaderInfoLog(vertexshader,300,0,infolog);
        cout<<"Vertex shader compilation failed\n";
    }
}
void OpenGL::fshader(){
    fragmentshader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentshader, 1,&fragment_shader,0);
    glCompileShader(fragmentshader);
    glGetShaderiv(fragmentshader, GL_COMPILE_STATUS, &s_Compile);
    if(!s_Compile){glGetShaderInfoLog(fragmentshader, 300, 0, infolog);
        cout<<"Fragment shader compilation failed"<<infolog<<'\n';
    }
}
void OpenGL::linkShaders(){
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram,  vertexshader);
    glAttachShader(shaderProgram, fragmentshader);
    glLinkProgram(shaderProgram);
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &s_Compile);
    if(!s_Compile){
        glGetProgramInfoLog(shaderProgram, 300, 0, infolog);
        cout<<"Linking shaders failed\n";
    }

    glUseProgram(shaderProgram);
    glDeleteShader(vertexshader);
    glDeleteShader(fragmentshader);
}
void OpenGL::genBuffer(){
    glGenBuffers(1, &Node_buffer);
    glGenVertexArrays(1,&VAO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER,Node_buffer);

    glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);
}
void OpenGL::delBuffer(){
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glDeleteVertexArrays(1,&VAO);
    glDeleteBuffers(1,&Node_buffer);
}
void OpenGL::render(std::vector<GLfloat>Vertices){

    genBuffer();

    glBufferData(GL_ARRAY_BUFFER,Vertices.size()*sizeof(GLfloat),Vertices.data(), GL_STATIC_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, 3);

    delBuffer();
}

void key_callback(GLFWwindow*window, int key, int scancode, int action, int mode){
    if((key == GLFW_KEY_ESCAPE||key == GLFW_KEY_SPACE) && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
}
float dot(vec2 v, vec2 w){
    float result = v.x() * w.x() + v.y() * w.y();
    return result;
}


const char *OpenGL::vertex_shader ="\
#version 330 core\n\
layout (location = 0) in vec3 position;\
out vec3 fragmentColor;\
void main(){\
gl_Position = vec4(position.x, position.y + 0.5, 0.0, 1.0);\
vec3 fragmentColor = vec3(0.0,0.0,0.0);\
}\
";

const char *OpenGL::fragment_shader = "\
#version 330 core\n\
in vec3 fragmentColor;\
out vec3 color;\
void main(){\
color = fragmentColor;\
}\
";
