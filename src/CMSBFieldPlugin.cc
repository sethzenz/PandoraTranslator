/**
 *  @file   LCContent/src/LCPlugins/LCBFieldPlugin.cc
 * 
 *  @brief  Implementation of the lc bfield plugin class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "PFCal/runPandora/interface/CMSBFieldCalculator.h"

using namespace pandora;

namespace cms_content
{

CMSBFieldCalculator::CMSBFieldCalculator(const float innerBField) :
    m_innerBField(innerBField)   
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CMSBFieldCalculator::GetBField(const CartesianVector &positionVector) const
{
    return m_innerBField;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CMSBFieldCalculator::Initialize()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CMSBFieldCalculator::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace cms_content
