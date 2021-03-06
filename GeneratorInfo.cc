#include "GeneratorInfo.h"

#include <algorithm>
#include <stdexcept>


pec::GeneratorInfo::GeneratorInfo() noexcept:
    processId(0),
    nominalWeight(0),
    pdfX(),  // the array is zeroed according to the C++03 standard
    pdfId(0),
    pdfQScale(0)
{}


void pec::GeneratorInfo::Reset()
{
    processId = 0;
    nominalWeight = 0;
    altWeights.clear();
    pdfId = 0;
    pdfX[0] = pdfX[1] = 0;
    pdfQScale = 0;
}


void pec::GeneratorInfo::SetProcessId(int processId_)
{
    processId = processId_;
}


void pec::GeneratorInfo::SetNominalWeight(float weight)
{
    nominalWeight = weight;
}


void pec::GeneratorInfo::AddAltWeight(float weight)
{
    altWeights.emplace_back(weight);
}


void pec::GeneratorInfo::SetPdfX(unsigned index, float x)
{
    // Check the index
    if (index > 1)
        throw std::logic_error("GeneratorInfo::SetPdfX: Illegal parton index.");
    
    
    // Check the desired fraction
    if (x < 0. or x > 1.)
        throw std::logic_error("GeneratorInfo::SetPdfX: The fraction must be in the range "
         "[0., 1.].");
    
    
    // Set the fraction
    pdfX[index] = x;
}


void pec::GeneratorInfo::SetPdfXs(float x1, float x2)
{
    SetPdfX(0, x1);
    SetPdfX(1, x2);
}


void pec::GeneratorInfo::SetPdfId(unsigned index, int id)
{
    // Check the index
    if (index > 1)
        throw std::logic_error("GeneratorInfo::SetPdfId: Illegal parton index.");
    
    
    // If gluons are encoded with their PDG ID code (as in Pythia 8), change it to zero
    if (id == 21)
        id = 0;
    
    
    // Check the parton ID code
    if (abs(id) > 5)
        throw std::logic_error("GeneratorInfo::SetPdfId: Illegal parton ID.");
    
    
    // Set ID of the specified parton. The other one is not changed. The given ID is increased by 5
    //to make it non-negative
    if (index == 0)
        pdfId = (pdfId - pdfId % 16) + UChar_t(id + 5);
    else
        pdfId = 16 * UChar_t(id + 5) + pdfId % 16;
}


void pec::GeneratorInfo::SetPdfIds(int id1, int id2)
{
    SetPdfId(0, id1);
    SetPdfId(1, id2);
}


void pec::GeneratorInfo::SetPdfQScale(float scale)
{
    pdfQScale = scale;
}


int pec::GeneratorInfo::ProcessId() const
{
    return processId;
}


float pec::GeneratorInfo::NominalWeight() const
{
    return nominalWeight;
}


std::vector<Float_t> const &pec::GeneratorInfo::AltWeights() const
{
    return altWeights;
}


float pec::GeneratorInfo::PdfX(unsigned index) const
{
    if (index > 1)
        throw std::logic_error("GeneratorInfo::PdfX: Illegal parton index.");
    
    return pdfX[index];
}


int pec::GeneratorInfo::PdfId(unsigned index) const
{
    if (index > 1)
        throw std::logic_error("GeneratorInfo::PdfId: Illegal parton index.");
    
    
    // Decode the parton ID
    int id;
    
    if (index == 0)
        id = int(pdfId) % 16 - 5;
        //^ Note the type conversion which is needed as otherwise the final result would be unsigned
    else
        id = int(pdfId) / 16 - 5;
    
    
    // Internally, gluons are encoded with code 0; return the PDG ID code for them instead
    if (id == 0)
        return 21;
    else
        return id;
}


float pec::GeneratorInfo::PdfQScale() const
{
    return pdfQScale;
}
