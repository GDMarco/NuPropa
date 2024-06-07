/* name of the plugin: myPlugin*/
%module(directors="1", threads="1", allprotected="1") nupropa

/* Exceptions required */
%include "attribute.i"
%include "exception.i"
%include "stdint.i"
%include "std_array.i"
%include "std_complex.i"
%include "std_container.i"
%include "std_iostream.i"
%include "std_list.i"
%include "std_map.i"
%include "std_set.i"
%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_vector.i"
%include "stl.i"
%include "typemaps.i"

/* Ignore list */
%ignore operator nupropa::Channels*;
%ignore operator nupropa::NeutrinoField*;

/*  define headers to include into the wrapper. These are the plugin headers
 *  and the CRPRopa headers.
 */
%{
#include "CRPropa.h"
#include "nupropa/NeutrinoNeutrinoInteraction.h"
#include "nupropa/NeutrinoAntineutrinoInteraction.h"
#include "nupropa/NeutrinoPhotonInteraction.h"
#include "nupropa/Channels.h"
#include "nupropa/NeutrinoBackground.h"

using namespace nupropa;
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/*Template for ref_ptr*/
%implicitconv crpropa::ref_ptr<nupropa::Channels>;
%template(ChannelsRefPtr) crpropa::ref_ptr<nupropa::Channels>;
%feature("director") nupropa::Channels;

/*Template for ref_ptr*/
%implicitconv crpropa::ref_ptr<nupropa::NeutrinoField>;
%template(NeutrinoFieldRefPtr) crpropa::ref_ptr<nupropa::NeutrinoField>;
%feature("director") nupropa::NeutrinoField;

/* include plugin parts to generate wrappers for */
%include "nupropa/NeutrinoNeutrinoInteraction.h"
%include "nupropa/NeutrinoAntineutrinoInteraction.h"
%include "nupropa/NeutrinoPhotonInteraction.h"
%include "nupropa/Channels.h"
%include "nupropa/NeutrinoBackground.h"




