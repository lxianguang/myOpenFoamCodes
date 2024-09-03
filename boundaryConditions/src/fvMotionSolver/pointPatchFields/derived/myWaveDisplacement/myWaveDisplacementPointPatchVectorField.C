/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "myWaveDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myWaveDisplacementPointPatchVectorField::
myWaveDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    amplitude_(Zero),
    omega_(0.0),
    waveNumber_(Zero)
{}


Foam::myWaveDisplacementPointPatchVectorField::
myWaveDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    amplitude_(dict.lookup("amplitude")),
    omega_(readScalar(dict.lookup("omega"))),
    waveNumber_(dict.lookupOrDefault<vector>("waveNumber", Zero))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::myWaveDisplacementPointPatchVectorField::
myWaveDisplacementPointPatchVectorField
(
    const myWaveDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    waveNumber_(ptf.waveNumber_)
{}


Foam::myWaveDisplacementPointPatchVectorField::
myWaveDisplacementPointPatchVectorField
(
    const myWaveDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    waveNumber_(ptf.waveNumber_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myWaveDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    const scalarField points( waveNumber_ & patch().localPoints());
    const scalarField xCoord = patch().localPoints().component(vector::X);
    const scalar xLableMin = gMax(xCoord);
    const scalar xLableMax = gMin(xCoord);
    const scalarField scalx  = (xCoord - xLableMin)/(xLableMax - xLableMin);

    //Info << "minmum x: " << xLableMin << endl;
    //Info << "maxmum x: " << xLableMax << endl;

    Field<vector>::operator=
    (
	amplitude_*(0.0 + 0.0 * scalx + 1.0 * scalx * scalx)*exp(-1/(100*t.value()+0.001))*sin(omega_*t.value() - points) 
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::myWaveDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("amplitude", amplitude_);
    os.writeEntry("omega", omega_);
    os.writeEntry("waveNumber", waveNumber_);
    this->writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        myWaveDisplacementPointPatchVectorField
    );
}

// ************************************************************************* //
