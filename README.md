**Setting Up Raylib and Compiling My Project**

### **1. Downloaded and Built Raylib**
- Cloned/downloaded Raylib source code.
- Navigated to the Raylib source directory:
  ```sh
  cd C:\raylib\raylib\src
  ```
- Compiled Raylib as a static library:
  ```sh
  mingw32-make PLATFORM=PLATFORM_DESKTOP
  ```
- This generated `libraylib.a` in `C:\raylib\raylib\src`.

### **2. Compiling My Project with Raylib**
- Navigated to my project directory:
  ```sh
  cd C:\Nash\Projects Folder\algorithmRunner
  ```
- Compiled using `g++` with static linking:
  ```sh
  g++ main.cpp -o game.exe -I "C:\raylib\raylib\src" -L "C:\raylib\raylib\src" -lraylib -lopengl32 -lgdi32 -lwinmm -static -g
  ./game.exe
  ```
- Generate static library
  ```sh
  cd C:\raylib\raylib\src
  mingw32-make PLATFORM=PLATFORM_DESKTOP
  ```

### **3. Debugging Issues**
- Initially got an error:
  ```
  could not convert 'InitWindow(...)' from 'void' to 'bool'
  ```
  **Solution:** `InitWindow()` returns void, so don’t use it inside `if (!InitWindow(...))`.
  
- Then encountered a crash with exit code `0xc0000139` (DLL entry point error).
  **Solution:** Ensured the game was linked against the freshly built static `libraylib.a`.
  
- Cleaned and rebuilt everything to make sure old objects weren’t interfering:
  ```sh
  mingw32-make clean
  mingw32-make PLATFORM=PLATFORM_DESKTOP
  g++ main.cpp -o game.exe -I "C:\raylib\raylib\src" -L "C:\raylib\raylib\src" -lraylib -lopengl32 -lgdi32 -lwinmm -static -g
  ```

### **4. Running the Game**
- Executed the compiled program:
  ```sh
  ./game.exe
  ```
- If any issues persist, run it with `gdb` for debugging:
  ```sh
  gdb game.exe
  run
  ```

### **5. Next Steps**
- If Raylib is updated, recompile it:
  ```sh
  mingw32-make clean
  mingw32-make PLATFORM=PLATFORM_DESKTOP
  ```
- Always use `-static` when linking Raylib to avoid missing DLL issues.
- If crashes persist, check missing dependencies using `Dependency Walker`.

**Now everything is working fine!** 🚀