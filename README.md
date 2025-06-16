# Implementación de BDF4 en OpenFOAM

Este repositorio está dedicado a la implementación del esquema de discretización temporal BDF4 en OpenFOAM, en el marco de la materia "Introducción a la Programación en OpenFOAM (2025)".

## Estructura del repositorio

```bash
BDF4_DdtScheme_OpenFOAM/
├── applications/
│   ├── solvers/
│   │   └── cosODEFoam/
│   │       ├── Make/
│   │       │   ├── files
│   │       │   └── options
│   │       ├── cosODEFoam.C
│   │       └── createFields.H
├── run/  # Casos para validar el orden del esquema. 
│   ├── cosODE/  # Casos 0D, campo uniforme en malla de 1 celda dependiente del tiempo
│   │   ├── case_bdf4_dt_.../
│   │   ├── BDF4_scheme_python_implementation.py
│   │   ├── generatePhiVsTimeCSVCase.py
│   │   └── postProcess.py
│   └── quickTest/  # Casos de advección pura con esquema QUICK + BDF4
│       ├── case_bdf4_quick_dt_..._dx_.../
│       ├── utility.py
│       └── postProcess.py
├── src/
│   └── finiteVolume/
│       ├── Make/
│       │   ├── files
│       │   └── options
│       └── finiteVolume/
│           └── ddtSchemes/
│               └── bdf4DdtScheme/
│                   ├── bdf4DdtScheme.C
│                   └── bdf4DdtScheme.H
├── docs/
│   └── INFORME_implementacion_esquema_bdf4_openfoam_FAZZARI.pdf # Memoria técnica de implementación y resultados
└── README.md
```
## Configuración y uso

1. **Clonar el repositorio**:

```bash
git clone https://github.com/rosen14/BDF4_DdtScheme_OpenFOAM.git
```

2. **Compilar el esquema**:

```bash
cd src/finiteVolume
wclean
wmake
```

3. Selecciona el esquema en tu fvSchemes
```bash
ddtSchemes
{
    default         bdf4;
}
```

4. Agrega la siguiente linea en tu controlDict:

```bash
libs (libmyBdf4DdtScheme);
```

## Contacto
Para preguntas o sugerencias sobre el proyecto, 
puedes contactar al mantenedor del repositorio a través de 
la página del proyecto en GitHub: [BDF4_DdtScheme_OpenFOAM](https://github.com/rosen14/BDF4_DdtScheme_OpenFOAM) o su email: rosendo.fazzari@ib.edu.ar

