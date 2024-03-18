{
    files = {
        "src/assignment-4.cpp"
    },
    values = {
        "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang",
        {
            "-Qunused-arguments",
            "-target",
            "arm64-apple-macos13.2",
            "-isysroot",
            "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX13.3.sdk",
            "-fvisibility=hidden",
            "-fvisibility-inlines-hidden",
            "-Wall",
            "-O3",
            "-std=c++20",
            "-Isrc/muni",
            "-isystem",
            "/opt/homebrew/Cellar/libomp/18.1.1/include",
            "-isystem",
            "/Users/rohankumar/.xmake/packages/s/spdlog/v1.13.0/ecca98d574b64b4ca0e391f301c3b293/include",
            "-isystem",
            "/Users/rohankumar/.xmake/packages/s/stb/2023.01.30/baa2aa3b50fd4294ba547a6c032c19ef/include",
            "-isystem",
            "/Users/rohankumar/.xmake/packages/s/stb/2023.01.30/baa2aa3b50fd4294ba547a6c032c19ef/include/stb",
            "-isystem",
            "/Users/rohankumar/.xmake/packages/l/linalg/v2.2/e58d89a1e9d448ca97233edda2f54d60/include",
            "-isystem",
            "/Users/rohankumar/.xmake/packages/t/tinyobjloader/v2.0.0rc13/138a7615eaf84e4f8c421e1d6ebdb490/include",
            "-Xpreprocessor",
            "-fopenmp",
            "-DNDEBUG"
        }
    },
    depfiles_gcc = "build/.objs/assignment-4/macosx/arm64/release/src/__cpp_assignment-4.cpp.cpp:   src/assignment-4.cpp src/muni/material.h src/muni/common.h   src/muni/math_helpers.h src/muni/camera.h src/muni/image.h   src/muni/obj_loader.h src/muni/triangle.h src/muni/ray_tracer.h   src/muni/sampler.h src/muni/scenes/box.h\
"
}