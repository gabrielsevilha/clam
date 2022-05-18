# Clam
## C Linear Algebra Math Library

  * C <br>
  * Unique Header <br>
  * Open Source <br>

```C

    Vector3 v = {0.0, 0.0, 0.0};
    
    Matrix4x4 p = createMatrix4x4(1.0);
    perspectiveMatrix(70.0f, width/height, 0.1f, 1000.0f, p);
    
    Matrix4x4 v = createMatrix4x4(1.0);
    lookatMatrix((Vector3){0.0,0.0,1.0},(Vector3){0.0,0.0,0.0},(Vector3){0.0,1.0,0.0}, p);
    
    //Send to opengl
    glUniformMatrix4fv(glGetUniformLocation(shader_program,"projection"), 1, GL_FALSE, p);

```
