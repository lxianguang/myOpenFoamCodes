/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 AUTHOR,AFFILIATION
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
    Foam::myFixedGradientFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition supplies a fixed gradient condition according to 
    the normal gradient, such that the patch values are calculated using:

        \f[
            x_p = x_c + \frac{\direction & \normal}{\Delta}
        \f]

    where
    \vartable
        x_p        | patch values
        x_c        | internal field values
        \direction | directional vector (user-specified)
        \normal    | normal gradient
        \Delta     | inverse distance from patch face centre to cell centre
    \endtable

Usage
    \table
        Property     | Description             | Required    | Default value
        direction    | direction               | yes         | (1,0,0)
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            myFixedGradient;
        direction       (1,0,0);
    }
    \endverbatim

SourceFiles
    myFixedGradientFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef myFixedGradientFvPatchScalarField_H
#define myFixedGradientFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class myFixedGradientFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class myFixedGradientFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data
        vector direction_;

public:

    //- Runtime type information
    TypeName("myFixedGradient");


    // Constructors

        //- Construct from patch and internal field
        myFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myFixedGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        myFixedGradientFvPatchScalarField
        (
            const myFixedGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        myFixedGradientFvPatchScalarField
        (
            const myFixedGradientFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new myFixedGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myFixedGradientFvPatchScalarField
        (
            const myFixedGradientFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new myFixedGradientFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            /*virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );*/


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
