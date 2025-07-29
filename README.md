This is a fluid simulation program. The smoke is simulated using Euler method (grid-based method), and rendered with volume ray casting. See [this video](https://www.youtube.com/watch?v=7UAFDXSXJu4&feature=youtu.be).


1. Prerequisites:
	
		opengl
		glew
		glfw
		Eigen

On Windows:
Step 1: Install dependencies via MSYS2
In MSYS2 MinGW 64-bit shell run:

		pacman -Syu
		pacman -S mingw-w64-x86_64-toolchain
		pacman -S mingw-w64-x86_64-glfw mingw-w64-x86_64-glew mingw-w64-x86_64-eigen3
		pacman -S mingw-w64-x86_64-opengl

Run to check shell: should print /mingw64/bin/g++

		which g++

2. How to run:
		
		1. make
		2. ./main

3. Control:

		Mouse:
		
			1. Change angle of view with mouse left key, and zoom with middle key.
			2. Select the Light and drag to change light position.

		Keyboard:
		
			1. R to reset the scene.
			2. S to switch between rendering and none rendering mode.
			3. W to toggle slices outline on/off.
			4. ESC to quit.
