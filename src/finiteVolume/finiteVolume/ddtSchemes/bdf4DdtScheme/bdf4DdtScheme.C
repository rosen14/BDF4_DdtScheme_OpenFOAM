/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2017,2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "bdf4DdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
scalar bdf4DdtScheme<Type>::deltaT_() const
{
    return mesh().time().deltaTValue();
}


template<class Type>
scalar bdf4DdtScheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0Value();
}


template<class Type>
template<class GeoField>
scalar bdf4DdtScheme<Type>::deltaT0_(const GeoField& vf) const
{
    if (mesh().time().timeIndex() < 2)
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}

template<class Type>
template<class GeoField>
scalar bdf4DdtScheme<Type>::deltaT00_(const GeoField& vf) const
{
    if (mesh().time().timeIndex() < 3)
    {
        return GREAT;
    }
    else
    {
        const scalar deltaT00 = deltaT0_();
        return deltaT00;
    }
}

template<class Type>
template<class GeoField>
scalar bdf4DdtScheme<Type>::deltaT000_(const GeoField& vf) const
{
    if (mesh().time().timeIndex() < 4)
    {
        return GREAT;
    }
    else
    {
        const scalar deltaT000 = deltaT0_();
        return deltaT000;
    }
}

template<class Type>
scalar bdf4DdtScheme<Type>::alpha_() const
{
    return deltaT_();
}

template<class Type>
template<class GeoField>
scalar bdf4DdtScheme<Type>::beta_(const GeoField& vf) const
{
    return deltaT_() + deltaT0_(vf);
}

template<class Type>
template<class GeoField>
scalar bdf4DdtScheme<Type>::gamma_(const GeoField& vf) const
{
    return deltaT_() + deltaT0_(vf) + deltaT00_(vf);
}

template<class Type>
template<class GeoField>
scalar bdf4DdtScheme<Type>::delta_(const GeoField& vf) const
{
    return deltaT_() + deltaT0_(vf) + deltaT00_(vf) + deltaT000_(vf);
}

template<class Type>
scalar bdf4DdtScheme<Type>::c0000_
(
    const scalar alphaC,
    const scalar beta,
    const scalar gamma,
    const scalar delta
) const
{
    return alphaC*beta*gamma /
           (delta*(delta - alphaC)*(delta - beta)*(delta - gamma));
}

template<class Type>
scalar bdf4DdtScheme<Type>::c000_
(
    const scalar alphaC,
    const scalar beta,
    const scalar gamma,
    const scalar delta
) const
{
    return alphaC * beta / 
           ((alphaC - gamma)*(gamma - beta)) * 
           (1/(delta - gamma) + 1/gamma);
}

template<class Type>
scalar bdf4DdtScheme<Type>::c00_
(
    const scalar alphaC,
    const scalar beta,
    const scalar gamma,
    const scalar delta
) const
{
    return alphaC / (beta - alphaC) * 
           (1/beta + 1/(gamma - beta)) + 
           alphaC * gamma / ((beta - alphaC)*(delta - gamma)) * 
           (1/(gamma - beta) - 1/(delta - beta));
}

template<class Type>
scalar bdf4DdtScheme<Type>::c0_
(
    const scalar alphaC,
    const scalar beta,
    const scalar gamma,
    const scalar delta
) const
{
    return beta * gamma / (delta - gamma) * 
               (1/((gamma - beta)*(gamma - alphaC)) + 
                1/((beta - alphaC)*(delta - beta))  - 
                1/((delta - alphaC)*(delta - beta)) - 
                1/((beta - alphaC)*(gamma - beta))) + 
           beta / (gamma - beta) *
               (1/(gamma - alphaC) - 1/(beta - alphaC)) -
           1/alphaC - 1/(beta - alphaC);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdf4DdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
    }

    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        ddtIOobject,
        mesh(),
        dimensioned<Type>(dt.dimensions()/dimTime, Zero),
        fvPatchFieldBase::calculatedType()
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdf4DdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar alphaC = alpha_();
    scalar beta  = beta_(vf);
    scalar gamma = gamma_(vf);
    scalar delta = delta_(vf);

    scalar coefft0000   = c0000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft000   = c000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft00 = c00_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft0  = c0_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft  = -(coefft0000 + coefft000 + coefft00 + coefft0);

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
        return tmp<GeometricField<Type, fvPatchField, volMesh>>(); // Retorno un valor dummy para evitar warning en compilación
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*vf
                  + coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                  + coefft000*vf.oldTime().oldTime().oldTime()
                  + coefft0000*vf.oldTime().oldTime().oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdf4DdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar alphaC = alpha_();
    scalar beta  = beta_(vf);
    scalar gamma = gamma_(vf);
    scalar delta = delta_(vf);

    scalar coefft0000   = c0000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft000   = c000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft00 = c00_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft0  = c0_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft  = -(coefft0000 + coefft000 + coefft00 + coefft0);

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
        return tmp<GeometricField<Type, fvPatchField, volMesh>>(); // Retorno un valor dummy para evitar warning en compilación
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                    coefft*vf
                  + coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                  + coefft000*vf.oldTime().oldTime().oldTime()
                  + coefft0000*vf.oldTime().oldTime().oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdf4DdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar alphaC = alpha_();
    scalar beta  = beta_(vf);
    scalar gamma = gamma_(vf);
    scalar delta = delta_(vf);

    scalar coefft0000   = c0000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft000   = c000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft00 = c00_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft0  = c0_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft  = -(coefft0000 + coefft000 + coefft00 + coefft0);

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
        return tmp<GeometricField<Type, fvPatchField, volMesh>>(); // Retorno un valor dummy para evitar warning en compilación
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*rho*vf
                  + coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                  + coefft000*rho.oldTime().oldTime().oldTime()*vf.oldTime().oldTime().oldTime()
                  + coefft0000*rho.oldTime().oldTime().oldTime().oldTime()*vf.oldTime().oldTime().oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
bdf4DdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar alphaC = alpha_();
    scalar beta  = beta_(vf);
    scalar gamma = gamma_(vf);
    scalar delta = delta_(vf);

    scalar coefft0000   = c0000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft000   = c000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft00 = c00_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft0  = c0_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft  = -(coefft0000 + coefft000 + coefft00 + coefft0);

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
        return tmp<GeometricField<Type, fvPatchField, volMesh>>(); // Retorno un valor dummy para evitar warning en compilación
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*alpha*rho*vf
                  + coefft0*
                        alpha.oldTime()*
                        rho.oldTime()*
                        vf.oldTime()
                  + coefft00*
                        alpha.oldTime().oldTime()*
                        rho.oldTime().oldTime()*
                        vf.oldTime().oldTime()
                  + coefft000*
                        alpha.oldTime().oldTime().oldTime()*
                        rho.oldTime().oldTime().oldTime()*
                        vf.oldTime().oldTime().oldTime()
                  + coefft0000*
                        alpha.oldTime().oldTime().oldTime().oldTime()*
                        rho.oldTime().oldTime().oldTime().oldTime()*
                        vf.oldTime().oldTime().oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type>>
bdf4DdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar alphaC = alpha_();
    scalar beta  = beta_(vf);
    scalar gamma = gamma_(vf);
    scalar delta = delta_(vf);

    scalar coefft0000   = c0000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft000   = c000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft00 = c00_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft0  = c0_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft  = -(coefft0000 + coefft000 + coefft00 + coefft0);

    fvm.diag() = (coefft*rDeltaT)*mesh().V();

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
        return tmp<fvMatrix<Type>>(); // Retorno un valor dummy para evitar warning en compilación
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
          - coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
          - coefft000*vf.oldTime().oldTime().oldTime().primitiveField()
          - coefft0000*vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
        );
    }
    
    Info << "------ BDF4 Coefficients ------" << nl;
    Info << "alphaC = " << alphaC << nl;
    Info << "beta   = " << beta   << nl;
    Info << "gamma  = " << gamma  << nl;
    Info << "delta  = " << delta  << nl;
    Info << "c0000  = " << coefft0000 << nl;
    Info << "c000   = " << coefft000  << nl;
    Info << "c00    = " << coefft00   << nl;
    Info << "c0     = " << coefft0    << nl;
    Info << "c      = " << coefft     << nl;
    
    /*
    Info << "------ Fields ------" << nl;
    Info << "vf: " << vf.primitiveField() << nl;
    Info << "vf.oldTime(): " << vf.oldTime().primitiveField() << nl;
    Info << "vf.oldTime().oldTime(): " << vf.oldTime().oldTime().primitiveField() << nl;
    Info << "vf.oldTime().oldTime().oldTime() " << vf.oldTime().oldTime().oldTime().primitiveField() << nl;
    Info << "vf.oldTime().oldTime().oldTime().oldTime() " << vf.oldTime().oldTime().oldTime().oldTime().primitiveField() << nl;
    Info << "--------------------------------" << endl;
    
    
    Info << "------ DeltaTimes ------" << nl;
    Info << "deltaT_(): " << deltaT_() << nl;
    Info << "deltaT0_(): " << deltaT0_(vf) << nl;
    Info << "deltaT00_(): " << deltaT00_(vf) << nl;
    Info << "deltaT000_() " << deltaT000_(vf) << nl;
    Info << "--------------------------------" << endl;
    */
    return tfvm;
}

template<class Type>
tmp<fvMatrix<Type>>
bdf4DdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar alphaC = alpha_();
    scalar beta  = beta_(vf);
    scalar gamma = gamma_(vf);
    scalar delta = delta_(vf);

    scalar coefft0000   = c0000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft000   = c000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft00 = c00_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft0  = c0_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft  = -(coefft0000 + coefft000 + coefft00 + coefft0);

    fvm.diag() = (coefft*rDeltaT*rho.value())*mesh().V();

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
        return tmp<fvMatrix<Type>>(); // Retorno un valor dummy para evitar warning en compilación
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*rho.value()*
        (
          - coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
          - coefft000*vf.oldTime().oldTime().oldTime().primitiveField()
          - coefft0000*vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
bdf4DdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar alphaC = alpha_();
    scalar beta  = beta_(vf);
    scalar gamma = gamma_(vf);
    scalar delta = delta_(vf);

    scalar coefft0000   = c0000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft000   = c000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft00 = c00_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft0  = c0_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft  = -(coefft0000 + coefft000 + coefft00 + coefft0);

    fvm.diag() = (coefft*rDeltaT)*rho.primitiveField()*mesh().V();

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
        return tmp<fvMatrix<Type>>(); // Retorno un valor dummy para evitar warning en compilación
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
          - coefft0*
                rho.oldTime().primitiveField()*
                vf.oldTime().primitiveField()
          - coefft00*
                rho.oldTime().oldTime().primitiveField()*
                vf.oldTime().oldTime().primitiveField()
          - coefft000*
                rho.oldTime().oldTime().oldTime().primitiveField()*
                vf.oldTime().oldTime().oldTime().primitiveField()
          - coefft0000*
                rho.oldTime().oldTime().oldTime().oldTime().primitiveField()*
                vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
bdf4DdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar alphaC = alpha_();
    scalar beta  = beta_(vf);
    scalar gamma = gamma_(vf);
    scalar delta = delta_(vf);

    scalar coefft0000   = c0000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft000   = c000_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft00 = c00_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft0  = c0_(alphaC,beta,gamma,delta)*deltaT_();
    scalar coefft  = -(coefft0000 + coefft000 + coefft00 + coefft0);

    fvm.diag() =
        (coefft*rDeltaT)*alpha.primitiveField()*rho.primitiveField()*mesh().V();

    if (mesh().moving())
    {
        FatalErrorInFunction
            << "The BDF4 time scheme is only implemented for static meshes." << nl
            << "Please use a static mesh or implement mesh motion support." << nl
            << abort(FatalError);
        return tmp<fvMatrix<Type>>(); // Retorno un valor dummy para evitar warning en compilación
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
          - coefft0
           *alpha.oldTime().primitiveField()
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()

          - coefft00
           *alpha.oldTime().oldTime().primitiveField()
           *rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()

          - coefft000
           *alpha.oldTime().oldTime().oldTime().primitiveField()
           *rho.oldTime().oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().oldTime().primitiveField()

          - coefft0000
           *alpha.oldTime().oldTime().oldTime().oldTime().primitiveField()
           *rho.oldTime().oldTime().oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<typename bdf4DdtScheme<Type>::fluxFieldType>
bdf4DdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    FatalErrorInFunction
        << "The method 'fvcDdtUfCorr' is not implemented for the BDF4 time scheme." << nl
        << "This scheme currently supports only 'fvmDdt'& 'fvcDdt' for static meshes." << nl
        << "If you require 'fvcDdtUfCorr' support, please implement it or use a different ddt scheme."  << nl
        << abort(FatalError);
    return tmp<fluxFieldType>(); // Retorno un valor dummy para evitar warning en compilación
}


template<class Type>
tmp<typename bdf4DdtScheme<Type>::fluxFieldType>
bdf4DdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    FatalErrorInFunction
        << "The method 'fvcDdtPhiCorr' is not implemented for the BDF4 time scheme." << nl
        << "This scheme currently supports only 'fvmDdt'& 'fvcDdt' for static meshes." << nl
        << "If you require 'fvcDdtPhiCorr' support, please implement it or use a different ddt scheme."  << nl
        << abort(FatalError);
    return tmp<fluxFieldType>(); // Retorno un valor dummy para evitar warning en compilación
}


template<class Type>
tmp<typename bdf4DdtScheme<Type>::fluxFieldType>
bdf4DdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    FatalErrorInFunction
        << "The method 'fvcDdtUfCorr' is not implemented for the BDF4 time scheme." << nl
        << "This scheme currently supports only 'fvmDdt'& 'fvcDdt' for static meshes." << nl
        << "If you require 'fvcDdtUfCorr' support, please implement it or use a different ddt scheme."  << nl
        << abort(FatalError);
    return tmp<fluxFieldType>(); // Retorno un valor dummy para evitar warning en compilación
}


template<class Type>
tmp<typename bdf4DdtScheme<Type>::fluxFieldType>
bdf4DdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    FatalErrorInFunction
        << "The method 'fvcDdtPhiCorr' is not implemented for the BDF4 time scheme." << nl
        << "This scheme currently supports only 'fvmDdt'& 'fvcDdt' for static meshes." << nl
        << "If you require 'fvcDdtPhiCorr' support, please implement it or use a different ddt scheme."  << nl
        << abort(FatalError);
    return tmp<fluxFieldType>(); // Retorno un valor dummy para evitar warning en compilación
}


template<class Type>
tmp<surfaceScalarField> bdf4DdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    FatalErrorInFunction
        << "The method 'meshPhi' is not implemented for the BDF4 time scheme." << nl
        << "This scheme currently supports only 'fvmDdt'& 'fvcDdt' for static meshes." << nl
        << "If you require 'meshPhi' support, please implement it or use a different ddt scheme."  << nl
        << abort(FatalError);
    return tmp<surfaceScalarField>(); // Retorno un valor dummy para evitar warning en compilación
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

makeFvDdtScheme(bdf4DdtScheme)
// ************************************************************************* //
