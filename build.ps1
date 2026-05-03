param(
    [string]$Compiler = "gfortran",
    [switch]$Test,
    [switch]$Clean,
    [string]$BuildDir = "build"
)

# ---- LASSI Fortran build script ----
$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
Set-Location $ScriptDir

if ($Clean) {
    if (Test-Path $BuildDir) { Remove-Item -Recurse -Force $BuildDir }
    Write-Host "Cleaned $BuildDir"
    return
}

# create build dir
if (-not (Test-Path $BuildDir)) { New-Item -ItemType Directory -Path $BuildDir | Out-Null }

# compile order matters because of module dependencies
$srcs = @(
    "src/lassi_kinds.f90",
    "src/lassi_geom.f90",
    "src/lassi_eos.f90",
    "src/lassi_friction.f90",
    "src/lassi_objects.f90",
    "src/lassi_grid.f90",
    "src/lassi_riemann.f90",
    "src/lassi_thomas.f90",
    "src/lassi_listmgmt.f90",
    "src/lassi_voidwave.f90",
    "src/lassi_press_mom.f90",
    "src/lassi_diag.f90",
    "src/lassi_io.f90",
    "src/lassi_main.f90"
)

# common flags
$cflags_release_gfortran = @("-O3","-march=native","-fopenmp","-Wno-uninitialized")
$cflags_debug_gfortran   = @("-O0","-g","-fbacktrace","-fcheck=all","-Wall")
$cflags_release_ifx      = @("/O3","/Qopenmp")
$cflags_debug_ifx        = @("/Od","/debug","/check:all","/traceback")

# pick compiler
$is_intel = $Compiler -match "^(ifx|ifort)$"
if ($is_intel) {
    $cflags = $cflags_release_ifx
    $modflag = "/module:$BuildDir"
    $outflag = "/exe:"
    $obj_ext = "obj"
} else {
    $cflags = $cflags_release_gfortran
    $modflag = "-J$BuildDir"
    $outflag = "-o "
    $obj_ext = "o"
}

# compile each source to .o
$objs = @()
foreach ($s in $srcs) {
    $name = [IO.Path]::GetFileNameWithoutExtension($s)
    $obj  = Join-Path $BuildDir "$name.$obj_ext"
    $objs += $obj
    Write-Host "[CC] $s -> $obj"
    if ($is_intel) {
        & $Compiler @cflags $modflag /c $s "/object:$obj"
    } else {
        & $Compiler @cflags $modflag -c $s -o $obj
    }
    if ($LASTEXITCODE -ne 0) { throw "Compile failed: $s" }
}

# link executable
$exe = Join-Path $BuildDir "lassi.exe"
Write-Host "[LD] -> $exe"
if ($is_intel) {
    & $Compiler @cflags $modflag $objs "/exe:$exe"
} else {
    & $Compiler @cflags $modflag $objs -o $exe
}
if ($LASTEXITCODE -ne 0) { throw "Link failed" }

if ($Test) {
    # build the unit test
    $test_obj = Join-Path $BuildDir "test_riemann.$obj_ext"
    Write-Host "[CC] tests/test_riemann.f90 -> $test_obj"
    if ($is_intel) {
        & $Compiler @cflags $modflag /c "tests/test_riemann.f90" "/object:$test_obj"
    } else {
        & $Compiler @cflags $modflag -c "tests/test_riemann.f90" -o $test_obj
    }
    if ($LASTEXITCODE -ne 0) { throw "Test compile failed" }

    # link test exe (needs the lassi modules' .o files except lassi_main.o)
    $core_objs = $objs | Where-Object { $_ -notmatch "lassi_main" }
    $test_exe = Join-Path $BuildDir "test_riemann.exe"
    if ($is_intel) {
        & $Compiler @cflags $modflag $test_obj $core_objs "/exe:$test_exe"
    } else {
        & $Compiler @cflags $modflag $test_obj $core_objs -o $test_exe
    }
    if ($LASTEXITCODE -ne 0) { throw "Test link failed" }
    Write-Host "Running test_riemann ..."
    & $test_exe
    if ($LASTEXITCODE -ne 0) { throw "Tests failed" }
}

Write-Host "Build OK."
