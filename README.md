# Clam
## C Linear Algebra Math Library

  * C <br>
  * Unique Header <br>
  * Open Source <br>

#### Usage Example:
```C
    Vector3 v = {0.0, 0.0, 0.0};
    
    Matrix4x4 p = INIT_MATRIX_4X4;
    perspectiveMatrix(70.0f, width/height, 0.1f, 1000.0f, p);
    
    Matrix4x4 a = INIT_MATRIX_4X4;
    lookatMatrix((Vector3){0.0,0.0,1.0},(Vector3){0.0,0.0,0.0},(Vector3){0.0,1.0,0.0}, a);
    
    //Send to opengl
    glUniformMatrix4fv(glGetUniformLocation(shader_program,"projection"), 1, GL_FALSE, p);
```
#### Matrix Multiplication Order:
```C
    //(MVP = projection x view x model)
    Matrix4x4 MVP = INIT_MATRIX_4X4;
    multiplyMatrix4x4(projection, view, MVP);
    multiplyMatrix4x4(MVP, model, MVP);
```
#### Transform:
```C
    Matrix4x4 m = INIT_MATRIX_4X4;
    
    identityMatrix4x4(m);
    
    translateMatrix4x4(m, (Vector3){5.0,10.0,15.0});
    rotateMatrix4x4(m, degToRad(90.0f), (Vector3){0.0,1.0,0.0});
    scaleMatrix4x4(m, (Vector3){2.0,2.0,2.0});
```
#### Quaternions:
```C
    Quaternion from = INIT_QUATERNION;
    Quaternion dest = anglesToQuaternion((Vector3){75.0f, 60.0f, 180.0f});

    Quaternion result = slerpQuaternion(from, dest, 0.5f);

    Matrix4x4 matrix;
    quaternionToMatrix4x4(result, matrix);
```
