{
    files = {
        "build/.objs/assignment-4/macosx/arm64/release/src/assignment-4.cpp.o"
    },
    values = {
        "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang++",
        {
            "-target",
            "arm64-apple-macos13.2",
            "-isysroot",
            "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX13.3.sdk",
            "-stdlib=libc++",
            "-lz",
            "-L/opt/homebrew/Cellar/libomp/18.1.1/lib",
            "-L/Users/rohankumar/.xmake/packages/t/tinyobjloader/v2.0.0rc13/138a7615eaf84e4f8c421e1d6ebdb490/lib",
            "-Wl,-x",
            "-Wl,-dead_strip",
            "-lomp",
            "-ltinyobjloader"
        }
    }
}