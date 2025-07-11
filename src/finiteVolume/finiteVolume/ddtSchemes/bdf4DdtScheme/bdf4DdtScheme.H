/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Class
    Foam::fv::bdf4DdtScheme

Group
    grpFvDdtSchemes

Description
    Fourth order backward-differencing ddt using the current and
    four previous time-step values.

SourceFiles
    bdf4DdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef bdf4DdtScheme_H
#define bdf4DdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class bdf4DdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class bdf4DdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous-previous time-step
        scalar deltaT00_() const;

        //- Return the previous-previous-previous time-step
        scalar deltaT000_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField&) const;

        //- Return the previous-previous time-step or GREAT if the old timestep field
        //  wasn't available in which case BDF2 ddt is used
        template<class GeoField>
        scalar deltaT00_(const GeoField&) const;

        //- Return the previous-previous-previous time-step or GREAT if the old timestep field
        //  wasn't available in which case BDF3 ddt is used
        template<class GeoField>
        scalar deltaT000_(const GeoField&) const;

        //- Return the alpha coefficient
        scalar alpha_() const;
        //- Return the beta coefficient
        template<class GeoField>
        scalar beta_(const GeoField&) const;
        //- Return the gamma coefficient
        template<class GeoField>
        scalar gamma_(const GeoField&) const;
        //- Return the delta coefficient
        template<class GeoField>
        scalar delta_(const GeoField&) const;

        //- Return the c0000 coefficient
        scalar c0000_
        (
            const scalar alphaC,
            const scalar beta,
            const scalar gamma,
            const scalar delta
        ) const;
        //- Return the c000 coefficient
        scalar c000_
        (
            const scalar alphaC,
            const scalar beta,
            const scalar gamma,
            const scalar delta
        ) const;
        //- Return the c00 coefficient
        scalar c00_
        (
            const scalar alphaC,
            const scalar beta,
            const scalar gamma,
            const scalar delta
        ) const;
        //- Return the c0 coefficient
        scalar c0_
        (
            const scalar alphaC,
            const scalar beta,
            const scalar gamma,
            const scalar delta
        ) const;


        //- No copy construct
        bdf4DdtScheme(const bdf4DdtScheme&) = delete;

        //- No copy assignment
        void operator=(const bdf4DdtScheme&) = delete;


public:

    //- Runtime type information
    TypeName("bdf4");


    // Constructors

        //- Construct from mesh
        bdf4DdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        bdf4DdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is)
        {
            if (is.good() && !is.eof())
            {
                is >> this->ddtPhiCoeff_;
            }

            // Ensure the old-old-time cell volumes are available
            // for moving meshes
            if (mesh.moving())
            {
                mesh.V00();
            }
        }


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensioned<Type>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


template<>
tmp<surfaceScalarField> bdf4DdtScheme<scalar>::fvcDdtUfCorr
(
    const GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf
);

template<>
tmp<surfaceScalarField> bdf4DdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> bdf4DdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> bdf4DdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
#ifdef NoRepository
    #include "bdf4DdtScheme.C"
#endif
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
