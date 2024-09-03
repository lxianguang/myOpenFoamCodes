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

\*---------------------------------------------------------------------------*/

#include "myFixedGradientFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myFixedGradientFvPatchScalarField::
myFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    direction_(1,0,0)
{
}


Foam::myFixedGradientFvPatchScalarField::
myFixedGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    direction_(dict.get<vector>("direction"))
{

    Info << "Using the myFixedGradient boundary condition" << endl;

    fixedValueFvPatchScalarField::evaluate();

    /*
    //Initialise with the value entry if evaluation is not possible
    readValueEntry(dict, IOobjectOption::MUST_READ);
    */
}


Foam::myFixedGradientFvPatchScalarField::
myFixedGradientFvPatchScalarField
(
    const myFixedGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    direction_(ptf.direction_)
{}


Foam::myFixedGradientFvPatchScalarField::
myFixedGradientFvPatchScalarField
(
    const myFixedGradientFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    direction_(ptf.direction_)
{}


Foam::myFixedGradientFvPatchScalarField::
myFixedGradientFvPatchScalarField
(
    const myFixedGradientFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    direction_(ptf.direction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*void Foam::myFixedGradientFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    fieldData_.autoMap(m);
}


void Foam::myFixedGradientFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const myFixedGradientFvPatchScalarField& tiptf =
        refCast<const myFixedGradientFvPatchScalarField>(ptf);

    fieldData_.rmap(tiptf.fieldData_, addr);
}*/


void Foam::myFixedGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField ngradient = direction_ & this->patch().nf();
    //const scalarField ngradient = this->patch().Sf().component(direction_)/this->patch().magSf();

    fixedValueFvPatchScalarField::operator==
    (
        this->patchInternalField() + ngradient/this->patch().deltaCoeffs()
    );


    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::myFixedGradientFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("direction", direction_);
    fvPatchScalarField::writeValueEntry(os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        myFixedGradientFvPatchScalarField
    );
}

// ************************************************************************* //
