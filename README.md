# Clam
## C Linear Algebra Math Library

  * C <br>
  * Unique Header <br>
  * Open Source <br>

```C

    Vector3 v = {0.0, 0.0, 0.0};
    
    Matrix4x4 p = INIT_MATRIX_4X4;
    perspectiveMatrix(70.0f, width/height, 0.1f, 1000.0f, p);
    
    Matrix4x4 a = INIT_MATRIX_4X4;
    lookatMatrix((Vector3){0.0,0.0,1.0},(Vector3){0.0,0.0,0.0},(Vector3){0.0,1.0,0.0}, a);
    
    //Send to opengl
    glUniformMatrix4fv(glGetUniformLocation(shader_program,"projection"), 1, GL_FALSE, p);
```

```C

    //Matrix multiplication order (MVP = projection x view x model)
    Matrix4x4 MVP = INIT_MATRIX_4X4;
    multiplyMatrix4x4(projection, view, MVP);
    multiplyMatrix4x4(MVP, model, MVP);
```
