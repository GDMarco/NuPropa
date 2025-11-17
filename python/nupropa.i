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
%ignore operator nupropa::ChannelsBundle*;
%ignore operator nupropa::NeutrinoField*;
%ignore operator nupropa::NeutrinoMixing*;
%ignore operator nupropa::RelativisticInteraction*;

/*  define headers to include into the wrapper. These are the plugin headers
 *  and the CRPRopa headers. [ParticleData is not a referenced pointer, just a class, that is called!]
 */
%{
#include "CRPropa.h"
#include "nupropa/NeutrinoNeutrinoInteraction.h"
#include "nupropa/NeutrinoAntineutrinoInteraction.h"
#include "nupropa/NeutrinoPhotonInteraction.h"
#include "nupropa/Channels.h"
#include "nupropa/ChannelsBundle.h"
#include "nupropa/NeutrinoField.h"
#include "nupropa/RelativisticInteraction.h"
#include "nupropa/ParticleData.h" 
#include "nupropa/NeutrinoMixing.h"
#include "nupropa/NeutrinoOscillation.h"

using namespace nupropa;
%}

%template(StringVector) std::vector<std::string>;
 %template(IntVector) std::vector<int>;

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/*Template for ref_ptr*/
%implicitconv crpropa::ref_ptr<nupropa::Channels>;
%template(ChannelsRefPtr) crpropa::ref_ptr<nupropa::Channels>;
%feature("director") nupropa::Channels;

%implicitconv crpropa::ref_ptr<nupropa::ChannelsBundle>;
%template(ChannelsBundleRefPtr) crpropa::ref_ptr<nupropa::ChannelsBundle>;
%feature("director") nupropa::ChannelsBundle;

%implicitconv crpropa::ref_ptr<nupropa::NeutrinoField>;
%template(NeutrinoFieldRefPtr) crpropa::ref_ptr<nupropa::NeutrinoField>;
%feature("director") nupropa::NeutrinoField;

%implicitconv crpropa::ref_ptr<nupropa::NeutrinoMixing>;
%template(NeutrinoMixingRefPtr) crpropa::ref_ptr<nupropa::NeutrinoMixing>;
%feature("director") nupropa::NeutrinoMixing;

%implicitconv crpropa::ref_ptr<nupropa::RelativisticInteraction>;
%template(RelativisticInteractionRefPtr) crpropa::ref_ptr<nupropa::RelativisticInteraction>;
%feature("director") nupropa::RelativisticInteraction;

/* see if they cause issues somewhere, probably not formally correct. It is also limited to give all the 4/5 args */
%extend nupropa::NeutrinoPhotonInteraction {
    // Full signature
    NeutrinoPhotonInteraction(crpropa::PhotonField *photonField,
                              nupropa::NeutrinoMixing *mixing,
                              bool haveSecondaries,
                              double limit) {
        return new nupropa::NeutrinoPhotonInteraction(
            crpropa::ref_ptr<crpropa::PhotonField>(photonField),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing),
            haveSecondaries,
            limit);
    }

    // photonField + mixing
    NeutrinoPhotonInteraction(crpropa::PhotonField *photonField,
                              nupropa::NeutrinoMixing *mixing) {
        return new nupropa::NeutrinoPhotonInteraction(
            crpropa::ref_ptr<crpropa::PhotonField>(photonField),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing));
    }

    // photonField + mixing + haveSecondaries
    NeutrinoPhotonInteraction(crpropa::PhotonField *photonField,
                              nupropa::NeutrinoMixing *mixing,
                              bool haveSecondaries) {
        return new nupropa::NeutrinoPhotonInteraction(
            crpropa::ref_ptr<crpropa::PhotonField>(photonField),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing),
            haveSecondaries);
    }
}

%extend nupropa::NeutrinoNeutrinoInteraction {
    // Full 4-argument version
    NeutrinoNeutrinoInteraction(nupropa::NeutrinoField *neutrinoField,
                                nupropa::NeutrinoMixing *mixing,
                                bool haveSecondaries,
                                double limit) {
        return new nupropa::NeutrinoNeutrinoInteraction(
            nupropa::ref_ptr<nupropa::NeutrinoField>(neutrinoField),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing),
            haveSecondaries,
            limit);
    }

    // 2-argument version
    NeutrinoNeutrinoInteraction(nupropa::NeutrinoField *neutrinoField,
                                nupropa::NeutrinoMixing *mixing) {
        return new nupropa::NeutrinoNeutrinoInteraction(
            nupropa::ref_ptr<nupropa::NeutrinoField>(neutrinoField),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing));
    }

    // 3-argument version
    NeutrinoNeutrinoInteraction(nupropa::NeutrinoField *neutrinoField,
                                nupropa::NeutrinoMixing *mixing,
                                bool haveSecondaries) {
        return new nupropa::NeutrinoNeutrinoInteraction(
            nupropa::ref_ptr<nupropa::NeutrinoField>(neutrinoField),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing),
            haveSecondaries);
    }
}

%extend nupropa::NeutrinoAntineutrinoInteraction {
    // Full 5-argument version
    NeutrinoAntineutrinoInteraction(nupropa::NeutrinoField *neutrinoField,
                                    nupropa::Channels *channels,
                                    nupropa::NeutrinoMixing *mixing,
                                    bool haveSecondaries,
                                    double limit) {
        return new nupropa::NeutrinoAntineutrinoInteraction(
            nupropa::ref_ptr<nupropa::NeutrinoField>(neutrinoField),
            crpropa::ref_ptr<nupropa::Channels>(channels),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing),
            haveSecondaries,
            limit);
    }

    // 3-argument version
    NeutrinoAntineutrinoInteraction(nupropa::NeutrinoField *neutrinoField,
                                    nupropa::Channels *channels,
                                    nupropa::NeutrinoMixing *mixing) {
        return new nupropa::NeutrinoAntineutrinoInteraction(
            nupropa::ref_ptr<nupropa::NeutrinoField>(neutrinoField),
            crpropa::ref_ptr<nupropa::Channels>(channels),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing));
    }

    // 4-argument version
    NeutrinoAntineutrinoInteraction(nupropa::NeutrinoField *neutrinoField,
                                    nupropa::Channels *channels,
                                    nupropa::NeutrinoMixing *mixing,
                                    bool haveSecondaries) {
        return new nupropa::NeutrinoAntineutrinoInteraction(
            nupropa::ref_ptr<nupropa::NeutrinoField>(neutrinoField),
            crpropa::ref_ptr<nupropa::Channels>(channels),
            crpropa::ref_ptr<nupropa::NeutrinoMixing>(mixing),
            haveSecondaries);
    }
}

/* include plugin parts to generate wrappers for */
%include "nupropa/NeutrinoNeutrinoInteraction.h"
%include "nupropa/NeutrinoAntineutrinoInteraction.h"
%include "nupropa/NeutrinoPhotonInteraction.h"
%include "nupropa/Channels.h"
%include "nupropa/ChannelsBundle.h"
%include "nupropa/NeutrinoField.h"
%include "nupropa/RelativisticInteraction.h"
%include "nupropa/ParticleData.h"
%include "nupropa/NeutrinoMixing.h"
%include "nupropa/NeutrinoOscillation.h"
%include "nupropa/NeutrinoField.h"




