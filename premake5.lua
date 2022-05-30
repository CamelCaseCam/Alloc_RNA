workspace "Alloc_RNA"
	architecture "x64"
	configurations
	{
		"Debug",
		"Release"
	}


OutputDir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"



project "Alloc_RNA"
    location "Alloc_RNA"
    kind "ConsoleApp"
    language "C++"

    targetdir ("bin/" .. OutputDir .. "/%{prj.name}")
	objdir ("bin/" .. OutputDir .. "/%{prj.name}")

	files
	{
		"%{prj.name}/src/Alloc_RNA.cpp",
		"%{prj.name}/src/FindDNA.cpp",
		"%{prj.name}/src/main.cpp",
		"%{prj.name}/src/Alloc_RNA.h",
		"%{prj.name}/src/FindDNA.h"
	}

	filter "system:windows"
		cppdialect "C++20"
		staticruntime "Off"
		systemversion "latest"

    filter "configurations:Debug"
		defines "Alloc_RNA_DEBUG"
		symbols "On"

    filter "configurations:Release"
		defines "Alloc_RNA_RELEASE"
		optimize "On"

    architecture "x64"
	configurations
	{
		"Debug",
		"Release"
	}
