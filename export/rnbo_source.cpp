/*******************************************************************************************************************
Cycling '74 License for Max-Generated Code for Export
Copyright (c) 2022 Cycling '74
The code that Max generates automatically and that end users are capable of exporting and using, and any
  associated documentation files (the “Software”) is a work of authorship for which Cycling '74 is the author
  and owner for copyright purposes.  A license is hereby granted, free of charge, to any person obtaining a
  copy of the Software (“Licensee”) to use, copy, modify, merge, publish, and distribute copies of the Software,
  and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The Software is licensed to Licensee only for non-commercial use. Users who wish to make commercial use of the
  Software must contact the copyright owner to determine if a license for commercial use is available, and the
  terms and conditions for same, which may include fees or royalties. For commercial use, please send inquiries
  to licensing@cycling74.com.  The determination of whether a use is commercial use or non-commercial use is based
  upon the use, not the user. The Software may be used by individuals, institutions, governments, corporations, or
  other business whether for-profit or non-profit so long as the use itself is not a commercialization of the
  materials or a use that generates or is intended to generate income, revenue, sales or profit.
The above copyright notice and this license shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
  THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
  CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.

Please see https://support.cycling74.com/hc/en-us/articles/10730637742483-RNBO-Export-Licensing-FAQ for additional information
*******************************************************************************************************************/

#include "RNBO_Common.h"
#include "RNBO_AudioSignal.h"

namespace RNBO {


#define floor(x) ((long)(x))

#if defined(__GNUC__) || defined(__clang__)
    #define RNBO_RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RNBO_RESTRICT __restrict
#endif

#define FIXEDSIZEARRAYINIT(...) { }

class rnbomatic : public PatcherInterfaceImpl {
public:

rnbomatic()
{
}

~rnbomatic()
{
}

rnbomatic* getTopLevelPatcher() {
    return this;
}

void cancelClockEvents()
{
}

template <typename T> void listquicksort(T& arr, T& sortindices, Int l, Int h, bool ascending) {
    if (l < h) {
        Int p = (Int)(this->listpartition(arr, sortindices, l, h, ascending));
        this->listquicksort(arr, sortindices, l, p - 1, ascending);
        this->listquicksort(arr, sortindices, p + 1, h, ascending);
    }
}

template <typename T> Int listpartition(T& arr, T& sortindices, Int l, Int h, bool ascending) {
    number x = arr[(Index)h];
    Int i = (Int)(l - 1);

    for (Int j = (Int)(l); j <= h - 1; j++) {
        bool asc = (bool)((bool)(ascending) && arr[(Index)j] <= x);
        bool desc = (bool)((bool)(!(bool)(ascending)) && arr[(Index)j] >= x);

        if ((bool)(asc) || (bool)(desc)) {
            i++;
            this->listswapelements(arr, i, j);
            this->listswapelements(sortindices, i, j);
        }
    }

    i++;
    this->listswapelements(arr, i, h);
    this->listswapelements(sortindices, i, h);
    return i;
}

template <typename T> void listswapelements(T& arr, Int a, Int b) {
    auto tmp = arr[(Index)a];
    arr[(Index)a] = arr[(Index)b];
    arr[(Index)b] = tmp;
}

inline number linearinterp(number frac, number x, number y) {
    return x + (y - x) * frac;
}

number maximum(number x, number y) {
    return (x < y ? y : x);
}

inline number safediv(number num, number denom) {
    return (denom == 0.0 ? 0.0 : num / denom);
}

number samplerate() {
    return this->sr;
}

MillisecondTime currenttime() {
    return this->_currentTime;
}

number tempo() {
    return this->getTopLevelPatcher()->globaltransport_getTempo();
}

number mstobeats(number ms) {
    return ms * this->tempo() * 0.008 / (number)480;
}

MillisecondTime sampstoms(number samps) {
    return samps * 1000 / this->sr;
}

ParameterIndex getParameterIndexForID(ConstCharPointer paramid) const {
    if (!stringCompare(paramid, "freq")) {
        return 0;
    }

    return INVALID_INDEX;
}

Index getNumMidiInputPorts() const {
    return 0;
}

void processMidiEvent(MillisecondTime , int , ConstByteArray , Index ) {}

Index getNumMidiOutputPorts() const {
    return 0;
}

void process(
    SampleValue ** inputs,
    Index numInputs,
    SampleValue ** outputs,
    Index numOutputs,
    Index n
) {
    RNBO_UNUSED(numInputs);
    RNBO_UNUSED(inputs);
    this->vs = n;
    this->updateTime(this->getEngine()->getCurrentTime());
    SampleValue * out1 = (numOutputs >= 1 && outputs[0] ? outputs[0] : this->dummyBuffer);
    SampleValue * out2 = (numOutputs >= 2 && outputs[1] ? outputs[1] : this->dummyBuffer);

    this->cycle_tilde_01_perform(
        this->cycle_tilde_01_frequency,
        this->cycle_tilde_01_phase_offset,
        this->signals[0],
        this->dummyBuffer,
        n
    );

    this->kink_tilde_01_perform(this->signals[0], this->kink_tilde_01_slope, out2, n);
    this->signalforwarder_01_perform(out2, out1, n);
    this->stackprotect_perform(n);
    this->globaltransport_advance();
    this->audioProcessSampleCount += this->vs;
}

void prepareToProcess(number sampleRate, Index maxBlockSize, bool force) {
    if (this->maxvs < maxBlockSize || !this->didAllocateSignals) {
        Index i;

        for (i = 0; i < 1; i++) {
            this->signals[i] = resizeSignal(this->signals[i], this->maxvs, maxBlockSize);
        }

        this->globaltransport_tempo = resizeSignal(this->globaltransport_tempo, this->maxvs, maxBlockSize);
        this->globaltransport_state = resizeSignal(this->globaltransport_state, this->maxvs, maxBlockSize);
        this->zeroBuffer = resizeSignal(this->zeroBuffer, this->maxvs, maxBlockSize);
        this->dummyBuffer = resizeSignal(this->dummyBuffer, this->maxvs, maxBlockSize);
        this->didAllocateSignals = true;
    }

    const bool sampleRateChanged = sampleRate != this->sr;
    const bool maxvsChanged = maxBlockSize != this->maxvs;
    const bool forceDSPSetup = sampleRateChanged || maxvsChanged || force;

    if (sampleRateChanged || maxvsChanged) {
        this->vs = maxBlockSize;
        this->maxvs = maxBlockSize;
        this->sr = sampleRate;
        this->invsr = 1 / sampleRate;
    }

    this->cycle_tilde_01_dspsetup(forceDSPSetup);
    this->globaltransport_dspsetup(forceDSPSetup);

    if (sampleRateChanged)
        this->onSampleRateChanged(sampleRate);
}

void setProbingTarget(MessageTag id) {
    switch (id) {
    default:
        this->setProbingIndex(-1);
        break;
    }
}

void setProbingIndex(ProbingIndex ) {}

Index getProbingChannels(MessageTag outletId) const {
    RNBO_UNUSED(outletId);
    return 0;
}

DataRef* getDataRef(DataRefIndex index)  {
    switch (index) {
    case 0:
        return addressOf(this->RNBODefaultMtofLookupTable256);
        break;
    case 1:
        return addressOf(this->RNBODefaultSinus);
        break;
    default:
        return nullptr;
    }
}

DataRefIndex getNumDataRefs() const {
    return 2;
}

void fillRNBODefaultMtofLookupTable256(DataRef& ref) {
    Float64BufferRef buffer;
    buffer = new Float64Buffer(ref);
    number bufsize = buffer->getSize();

    for (Index i = 0; i < bufsize; i++) {
        number midivalue = -256. + (number)512. / (bufsize - 1) * i;
        buffer[i] = rnbo_exp(.057762265 * (midivalue - 69.0));
    }
}

void fillRNBODefaultSinus(DataRef& ref) {
    Float64BufferRef buffer;
    buffer = new Float64Buffer(ref);
    number bufsize = buffer->getSize();

    for (Index i = 0; i < bufsize; i++) {
        buffer[i] = rnbo_cos(i * 3.14159265358979323846 * 2. / bufsize);
    }
}

void fillDataRef(DataRefIndex index, DataRef& ref) {
    switch (index) {
    case 0:
        this->fillRNBODefaultMtofLookupTable256(ref);
        break;
    case 1:
        this->fillRNBODefaultSinus(ref);
        break;
    }
}

void processDataViewUpdate(DataRefIndex index, MillisecondTime time) {
    this->updateTime(time);

    if (index == 0) {
        this->mtof_01_innerMtoF_buffer = new Float64Buffer(this->RNBODefaultMtofLookupTable256);
    }

    if (index == 1) {
        this->cycle_tilde_01_buffer = new Float64Buffer(this->RNBODefaultSinus);
        this->cycle_tilde_01_bufferUpdated();
    }
}

void initialize() {
    this->RNBODefaultMtofLookupTable256 = initDataRef("RNBODefaultMtofLookupTable256", true, nullptr);
    this->RNBODefaultSinus = initDataRef("RNBODefaultSinus", true, nullptr);
    this->assign_defaults();
    this->setState();
    this->RNBODefaultMtofLookupTable256->setIndex(0);
    this->mtof_01_innerMtoF_buffer = new Float64Buffer(this->RNBODefaultMtofLookupTable256);
    this->RNBODefaultSinus->setIndex(1);
    this->cycle_tilde_01_buffer = new Float64Buffer(this->RNBODefaultSinus);
    this->initializeObjects();
    this->allocateDataRefs();
    this->startup();
}

Index getIsMuted()  {
    return this->isMuted;
}

void setIsMuted(Index v)  {
    this->isMuted = v;
}

Index getPatcherSerial() const {
    return 0;
}

void getState(PatcherStateInterface& ) {}

void setState() {}

void getPreset(PatcherStateInterface& preset) {
    preset["__presetid"] = "rnbo";
    this->param_01_getPresetValue(getSubState(preset, "freq"));
}

void setPreset(MillisecondTime time, PatcherStateInterface& preset) {
    this->updateTime(time);
    this->param_01_setPresetValue(getSubState(preset, "freq"));
}

void processTempoEvent(MillisecondTime time, Tempo tempo) {
    this->updateTime(time);

    if (this->globaltransport_setTempo(tempo, false))
        {}
}

void processTransportEvent(MillisecondTime time, TransportState state) {
    this->updateTime(time);

    if (this->globaltransport_setState(state, false))
        {}
}

void processBeatTimeEvent(MillisecondTime time, BeatTime beattime) {
    this->updateTime(time);

    if (this->globaltransport_setBeatTime(beattime, false))
        {}
}

void onSampleRateChanged(double ) {}

void processTimeSignatureEvent(MillisecondTime time, int numerator, int denominator) {
    this->updateTime(time);

    if (this->globaltransport_setTimeSignature(numerator, denominator, false))
        {}
}

void setParameterValue(ParameterIndex index, ParameterValue v, MillisecondTime time) {
    this->updateTime(time);

    switch (index) {
    case 0:
        this->param_01_value_set(v);
        break;
    }
}

void processParameterEvent(ParameterIndex index, ParameterValue value, MillisecondTime time) {
    this->setParameterValue(index, value, time);
}

void processNormalizedParameterEvent(ParameterIndex index, ParameterValue value, MillisecondTime time) {
    this->setParameterValueNormalized(index, value, time);
}

ParameterValue getParameterValue(ParameterIndex index)  {
    switch (index) {
    case 0:
        return this->param_01_value;
    default:
        return 0;
    }
}

ParameterIndex getNumSignalInParameters() const {
    return 0;
}

ParameterIndex getNumSignalOutParameters() const {
    return 0;
}

ParameterIndex getNumParameters() const {
    return 1;
}

ConstCharPointer getParameterName(ParameterIndex index) const {
    switch (index) {
    case 0:
        return "freq";
    default:
        return "bogus";
    }
}

ConstCharPointer getParameterId(ParameterIndex index) const {
    switch (index) {
    case 0:
        return "freq";
    default:
        return "bogus";
    }
}

void getParameterInfo(ParameterIndex index, ParameterInfo * info) const {
    {
        switch (index) {
        case 0:
            info->type = ParameterTypeNumber;
            info->initialValue = 220;
            info->min = 100;
            info->max = 2000;
            info->exponent = 1;
            info->steps = 0;
            info->debug = false;
            info->saveable = true;
            info->transmittable = true;
            info->initialized = true;
            info->visible = true;
            info->displayName = "";
            info->unit = "";
            info->ioType = IOTypeUndefined;
            info->signalIndex = INVALID_INDEX;
            break;
        }
    }
}

void sendParameter(ParameterIndex index, bool ignoreValue) {
    this->getEngine()->notifyParameterValueChanged(index, (ignoreValue ? 0 : this->getParameterValue(index)), ignoreValue);
}

ParameterValue applyStepsToNormalizedParameterValue(ParameterValue normalizedValue, int steps) const {
    if (steps == 1) {
        if (normalizedValue > 0) {
            normalizedValue = 1.;
        }
    } else {
        ParameterValue oneStep = (number)1. / (steps - 1);
        ParameterValue numberOfSteps = rnbo_fround(normalizedValue / oneStep * 1 / (number)1) * (number)1;
        normalizedValue = numberOfSteps * oneStep;
    }

    return normalizedValue;
}

ParameterValue convertToNormalizedParameterValue(ParameterIndex index, ParameterValue value) const {
    switch (index) {
    case 0:
        {
            value = (value < 100 ? 100 : (value > 2000 ? 2000 : value));
            ParameterValue normalizedValue = (value - 100) / (2000 - 100);
            return normalizedValue;
        }
    default:
        return value;
    }
}

ParameterValue convertFromNormalizedParameterValue(ParameterIndex index, ParameterValue value) const {
    value = (value < 0 ? 0 : (value > 1 ? 1 : value));

    switch (index) {
    case 0:
        {
            value = (value < 0 ? 0 : (value > 1 ? 1 : value));

            {
                return 100 + value * (2000 - 100);
            }
        }
    default:
        return value;
    }
}

ParameterValue constrainParameterValue(ParameterIndex index, ParameterValue value) const {
    switch (index) {
    case 0:
        return this->param_01_value_constrain(value);
    default:
        return value;
    }
}

void scheduleParamInit(ParameterIndex index, Index order) {
    this->paramInitIndices->push(index);
    this->paramInitOrder->push(order);
}

void processParamInitEvents() {
    this->listquicksort(
        this->paramInitOrder,
        this->paramInitIndices,
        0,
        (int)(this->paramInitOrder->length - 1),
        true
    );

    for (Index i = 0; i < this->paramInitOrder->length; i++) {
        this->getEngine()->scheduleParameterChange(
            this->paramInitIndices[i],
            this->getParameterValue(this->paramInitIndices[i]),
            0
        );
    }
}

void processClockEvent(MillisecondTime , ClockId , bool , ParameterValue ) {}

void processOutletAtCurrentTime(EngineLink* , OutletIndex , ParameterValue ) {}

void processOutletEvent(
    EngineLink* sender,
    OutletIndex index,
    ParameterValue value,
    MillisecondTime time
) {
    this->updateTime(time);
    this->processOutletAtCurrentTime(sender, index, value);
}

void processNumMessage(MessageTag , MessageTag , MillisecondTime , number ) {}

void processListMessage(MessageTag , MessageTag , MillisecondTime , const list& ) {}

void processBangMessage(MessageTag , MessageTag , MillisecondTime ) {}

MessageTagInfo resolveTag(MessageTag tag) const {
    switch (tag) {

    }

    return "";
}

MessageIndex getNumMessages() const {
    return 0;
}

const MessageInfo& getMessageInfo(MessageIndex index) const {
    switch (index) {

    }

    return NullMessageInfo;
}

protected:

void param_01_value_set(number v) {
    v = this->param_01_value_constrain(v);
    this->param_01_value = v;
    this->sendParameter(0, false);

    if (this->param_01_value != this->param_01_lastValue) {
        this->getEngine()->presetTouched();
        this->param_01_lastValue = this->param_01_value;
    }

    {
        list converted = {v};
        this->ftom_01_frequency_set(converted);
    }
}

number msToSamps(MillisecondTime ms, number sampleRate) {
    return ms * sampleRate * 0.001;
}

MillisecondTime sampsToMs(SampleIndex samps) {
    return samps * (this->invsr * 1000);
}

Index getMaxBlockSize() const {
    return this->maxvs;
}

number getSampleRate() const {
    return this->sr;
}

bool hasFixedVectorSize() const {
    return false;
}

Index getNumInputChannels() const {
    return 0;
}

Index getNumOutputChannels() const {
    return 2;
}

void allocateDataRefs() {
    this->mtof_01_innerMtoF_buffer->requestSize(65536, 1);
    this->mtof_01_innerMtoF_buffer->setSampleRate(this->sr);
    this->cycle_tilde_01_buffer->requestSize(16384, 1);
    this->cycle_tilde_01_buffer->setSampleRate(this->sr);
    this->mtof_01_innerMtoF_buffer = this->mtof_01_innerMtoF_buffer->allocateIfNeeded();

    if (this->RNBODefaultMtofLookupTable256->hasRequestedSize()) {
        if (this->RNBODefaultMtofLookupTable256->wantsFill())
            this->fillRNBODefaultMtofLookupTable256(this->RNBODefaultMtofLookupTable256);

        this->getEngine()->sendDataRefUpdated(0);
    }

    this->cycle_tilde_01_buffer = this->cycle_tilde_01_buffer->allocateIfNeeded();

    if (this->RNBODefaultSinus->hasRequestedSize()) {
        if (this->RNBODefaultSinus->wantsFill())
            this->fillRNBODefaultSinus(this->RNBODefaultSinus);

        this->getEngine()->sendDataRefUpdated(1);
    }
}

void initializeObjects() {
    this->mtof_01_innerScala_init();
    this->mtof_01_init();
    this->ftom_01_innerScala_init();
    this->ftom_01_init();
}

void sendOutlet(OutletIndex index, ParameterValue value) {
    this->getEngine()->sendOutlet(this, index, value);
}

void startup() {
    this->updateTime(this->getEngine()->getCurrentTime());

    {
        this->scheduleParamInit(0, 0);
    }

    this->processParamInitEvents();
}

static number param_01_value_constrain(number v) {
    v = (v > 2000 ? 2000 : (v < 100 ? 100 : v));
    return v;
}

void cycle_tilde_01_frequency_set(number v) {
    this->cycle_tilde_01_frequency = v;
}

void cycle_tilde_01_phase_offset_set(number v) {
    this->cycle_tilde_01_phase_offset = v;
}

void mtof_01_out_set(const list& v) {
    {
        if (v->length > 1)
            this->cycle_tilde_01_phase_offset_set(v[1]);

        number converted = (v->length > 0 ? v[0] : 0);
        this->cycle_tilde_01_frequency_set(converted);
    }
}

void mtof_01_midivalue_set(const list& v) {
    this->mtof_01_midivalue = jsCreateListCopy(v);
    list tmp = list();

    for (int i = 0; i < this->mtof_01_midivalue->length; i++) {
        tmp->push(
            this->mtof_01_innerMtoF_next(this->mtof_01_midivalue[(Index)i], this->mtof_01_base)
        );
    }

    this->mtof_01_out_set(tmp);
}

void ftom_01_out_set(const list& v) {
    this->mtof_01_midivalue_set(v);
}

void ftom_01_frequency_set(const list& v) {
    this->ftom_01_frequency = jsCreateListCopy(v);
    list tmp = list();

    for (int i = 0; i < this->ftom_01_frequency->length; i++) {
        number v = this->ftom_01_innerFtoM_next(this->ftom_01_frequency[(Index)i], this->ftom_01_base);

        if (1 == 0 || v != -999.0) {
            tmp->push(v);
        }
    }

    if (1 == 0 || tmp->length > 0) {
        this->ftom_01_out_set(tmp);
    }
}

void cycle_tilde_01_perform(
    number frequency,
    number phase_offset,
    Sample * out1,
    Sample * out2,
    Index n
) {
    RNBO_UNUSED(phase_offset);
    auto __cycle_tilde_01_f2i = this->cycle_tilde_01_f2i;
    auto __cycle_tilde_01_phasei = this->cycle_tilde_01_phasei;
    Index i;

    for (i = 0; i < n; i++) {
        {
            uint32_t uint_phase;

            {
                {
                    uint_phase = __cycle_tilde_01_phasei;
                }
            }

            uint32_t idx = (uint32_t)(uint32_rshift(uint_phase, 18));
            number frac = ((BinOpInt)((UBinOpInt)uint_phase & (UBinOpInt)262143)) * 3.81471181759574e-6;
            number y0 = this->cycle_tilde_01_buffer[(Index)idx];
            number y1 = this->cycle_tilde_01_buffer[(Index)((UBinOpInt)(idx + 1) & (UBinOpInt)16383)];
            number y = y0 + frac * (y1 - y0);

            {
                uint32_t pincr = (uint32_t)(uint32_trunc(frequency * __cycle_tilde_01_f2i));
                __cycle_tilde_01_phasei = uint32_add(__cycle_tilde_01_phasei, pincr);
            }

            out1[(Index)i] = y;
            out2[(Index)i] = uint_phase * 0.232830643653869629e-9;
            continue;
        }
    }

    this->cycle_tilde_01_phasei = __cycle_tilde_01_phasei;
}

void kink_tilde_01_perform(const Sample * x, number slope, Sample * out1, Index n) {
    RNBO_UNUSED(slope);
    Index i;

    for (i = 0; i < n; i++) {
        number clippedSlope = this->maximum(1e-9, 0.5);
        number cross = (number)(this->safediv(0.5, clippedSlope));

        if (x[(Index)i] > cross) {
            number s2 = (number)((number)0.5 / (1.0 - cross));
            number ic = (number)(1.0 - s2);
            out1[(Index)i] = x[(Index)i] * s2 + ic;
            continue;
        } else {
            out1[(Index)i] = x[(Index)i] * clippedSlope;
            continue;
        }
    }
}

void signalforwarder_01_perform(const Sample * input, Sample * output, Index n) {
    for (Index i = 0; i < n; i++) {
        output[(Index)i] = input[(Index)i];
    }
}

void stackprotect_perform(Index n) {
    RNBO_UNUSED(n);
    auto __stackprotect_count = this->stackprotect_count;
    __stackprotect_count = 0;
    this->stackprotect_count = __stackprotect_count;
}

number mtof_01_innerMtoF_next(number midivalue, number tuning) {
    if (midivalue == this->mtof_01_innerMtoF_lastInValue && tuning == this->mtof_01_innerMtoF_lastTuning)
        return this->mtof_01_innerMtoF_lastOutValue;

    this->mtof_01_innerMtoF_lastInValue = midivalue;
    this->mtof_01_innerMtoF_lastTuning = tuning;
    number result = 0;

    {
        result = rnbo_exp(.057762265 * (midivalue - 69.0));
    }

    this->mtof_01_innerMtoF_lastOutValue = tuning * result;
    return this->mtof_01_innerMtoF_lastOutValue;
}

void mtof_01_innerMtoF_reset() {
    this->mtof_01_innerMtoF_lastInValue = 0;
    this->mtof_01_innerMtoF_lastOutValue = 0;
    this->mtof_01_innerMtoF_lastTuning = 0;
}

void mtof_01_innerScala_mid(int v) {
    this->mtof_01_innerScala_kbmMid = v;
    this->mtof_01_innerScala_updateRefFreq();
}

void mtof_01_innerScala_ref(int v) {
    this->mtof_01_innerScala_kbmRefNum = v;
    this->mtof_01_innerScala_updateRefFreq();
}

void mtof_01_innerScala_base(number v) {
    this->mtof_01_innerScala_kbmRefFreq = v;
    this->mtof_01_innerScala_updateRefFreq();
}

void mtof_01_innerScala_init() {
    list sclValid = {
        12,
        100,
        0,
        200,
        0,
        300,
        0,
        400,
        0,
        500,
        0,
        600,
        0,
        700,
        0,
        800,
        0,
        900,
        0,
        1000,
        0,
        1100,
        0,
        2,
        1
    };

    this->mtof_01_innerScala_updateScale(sclValid);
}

void mtof_01_innerScala_update(list scale, list map) {
    if (scale->length > 0) {
        this->mtof_01_innerScala_updateScale(scale);
    }

    if (map->length > 0) {
        this->mtof_01_innerScala_updateMap(map);
    }
}

number mtof_01_innerScala_mtof(number note) {
    if ((bool)(this->mtof_01_innerScala_lastValid) && this->mtof_01_innerScala_lastNote == note) {
        return this->mtof_01_innerScala_lastFreq;
    }

    array<int, 2> degoct = this->mtof_01_innerScala_applyKBM(note);
    number out = 0;

    if (degoct[1] > 0) {
        out = this->mtof_01_innerScala_applySCL(degoct[0], fract(note), this->mtof_01_innerScala_refFreq);
    }

    this->mtof_01_innerScala_updateLast(note, out);
    return out;
}

number mtof_01_innerScala_ftom(number hz) {
    if (hz <= 0.0) {
        return 0.0;
    }

    if ((bool)(this->mtof_01_innerScala_lastValid) && this->mtof_01_innerScala_lastFreq == hz) {
        return this->mtof_01_innerScala_lastNote;
    }

    array<number, 2> df = this->mtof_01_innerScala_hztodeg(hz);
    int degree = (int)(df[0]);
    number frac = df[1];
    number out = 0;

    if (this->mtof_01_innerScala_kbmSize == 0) {
        out = this->mtof_01_innerScala_kbmMid + degree;
    } else {
        array<int, 2> octdeg = this->mtof_01_innerScala_octdegree(degree, this->mtof_01_innerScala_kbmOctaveDegree);
        number oct = (number)(octdeg[0]);
        int index = (int)(octdeg[1]);
        Index entry = 0;

        for (Index i = 0; i < this->mtof_01_innerScala_kbmMapSize; i++) {
            if (index == this->mtof_01_innerScala_kbmValid[(Index)(i + this->mtof_01_innerScala_KBM_MAP_OFFSET)]) {
                entry = i;
                break;
            }
        }

        out = oct * this->mtof_01_innerScala_kbmSize + entry + this->mtof_01_innerScala_kbmMid;
    }

    out = out + frac;
    this->mtof_01_innerScala_updateLast(out, hz);
    return this->mtof_01_innerScala_lastNote;
}

int mtof_01_innerScala_updateScale(list scl) {
    if (scl->length > 1 && scl[0] * 2 + 1 == scl->length) {
        this->mtof_01_innerScala_lastValid = false;
        this->mtof_01_innerScala_sclExpMul = {};
        number last = 1;

        for (Index i = 1; i < scl->length; i += 2) {
            const number c = (const number)(scl[(Index)(i + 0)]);
            const number d = (const number)(scl[(Index)(i + 1)]);

            if (d <= 0) {
                last = c / (number)1200;
            } else {
                last = rnbo_log2(c / d);
            }

            this->mtof_01_innerScala_sclExpMul->push(last);
        }

        this->mtof_01_innerScala_sclOctaveMul = last;
        this->mtof_01_innerScala_sclEntryCount = (int)(this->mtof_01_innerScala_sclExpMul->length);
        this->mtof_01_innerScala_updateRefFreq();
        return 1;
    }

    return 0;
}

int mtof_01_innerScala_updateMap(list kbm) {
    if (kbm->length == 1 && kbm[0] == 0.0) {
        kbm = {0.0, 0.0, 0.0, 60.0, 69.0, 440.0};
    }

    if (kbm->length >= 6 && kbm[0] >= 0.0) {
        this->mtof_01_innerScala_lastValid = false;
        Index size = (Index)(kbm[0]);
        int octave = 12;

        if (kbm->length > 6) {
            octave = (int)(kbm[6]);
        }

        if (size > 0 && kbm->length < this->mtof_01_innerScala_KBM_MAP_OFFSET) {
            return 0;
        }

        this->mtof_01_innerScala_kbmSize = (int)(size);
        this->mtof_01_innerScala_kbmMin = (int)(kbm[1]);
        this->mtof_01_innerScala_kbmMax = (int)(kbm[2]);
        this->mtof_01_innerScala_kbmMid = (int)(kbm[3]);
        this->mtof_01_innerScala_kbmRefNum = (int)(kbm[4]);
        this->mtof_01_innerScala_kbmRefFreq = kbm[5];
        this->mtof_01_innerScala_kbmOctaveDegree = octave;
        this->mtof_01_innerScala_kbmValid = kbm;
        this->mtof_01_innerScala_kbmMapSize = (kbm->length - this->mtof_01_innerScala_KBM_MAP_OFFSET > kbm->length ? kbm->length : (kbm->length - this->mtof_01_innerScala_KBM_MAP_OFFSET < 0 ? 0 : kbm->length - this->mtof_01_innerScala_KBM_MAP_OFFSET));
        this->mtof_01_innerScala_updateRefFreq();
        return 1;
    }

    return 0;
}

void mtof_01_innerScala_updateLast(number note, number freq) {
    this->mtof_01_innerScala_lastValid = true;
    this->mtof_01_innerScala_lastNote = note;
    this->mtof_01_innerScala_lastFreq = freq;
}

array<number, 2> mtof_01_innerScala_hztodeg(number hz) {
    number hza = rnbo_abs(hz);

    number octave = rnbo_floor(
        rnbo_log2(hza / this->mtof_01_innerScala_refFreq) / this->mtof_01_innerScala_sclOctaveMul
    );

    int i = 0;
    number frac = 0;
    number n = 0;

    for (; i < this->mtof_01_innerScala_sclEntryCount; i++) {
        number c = this->mtof_01_innerScala_applySCLOctIndex(octave, i + 0, 0.0, this->mtof_01_innerScala_refFreq);
        n = this->mtof_01_innerScala_applySCLOctIndex(octave, i + 1, 0.0, this->mtof_01_innerScala_refFreq);

        if (c <= hza && hza < n) {
            if (c != hza) {
                frac = rnbo_log2(hza / c) / rnbo_log2(n / c);
            }

            break;
        }
    }

    if (i == this->mtof_01_innerScala_sclEntryCount && n != hza) {
        number c = n;
        n = this->mtof_01_innerScala_applySCLOctIndex(octave + 1, 0, 0.0, this->mtof_01_innerScala_refFreq);
        frac = rnbo_log2(hza / c) / rnbo_log2(n / c);
    }

    number deg = i + octave * this->mtof_01_innerScala_sclEntryCount;

    {
        deg = rnbo_fround((deg + frac) * 1 / (number)1) * 1;
        frac = 0.0;
    }

    return {deg, frac};
}

array<int, 2> mtof_01_innerScala_octdegree(int degree, int count) {
    int octave = 0;
    int index = 0;

    if (degree < 0) {
        octave = -(1 + (-1 - degree) / count);
        index = -degree % count;

        if (index > 0) {
            index = count - index;
        }
    } else {
        octave = degree / count;
        index = degree % count;
    }

    return {octave, index};
}

array<int, 2> mtof_01_innerScala_applyKBM(number note) {
    if ((this->mtof_01_innerScala_kbmMin == this->mtof_01_innerScala_kbmMax && this->mtof_01_innerScala_kbmMax == 0) || (note >= this->mtof_01_innerScala_kbmMin && note <= this->mtof_01_innerScala_kbmMax)) {
        int degree = (int)(rnbo_floor(note - this->mtof_01_innerScala_kbmMid));

        if (this->mtof_01_innerScala_kbmSize == 0) {
            return {degree, 1};
        }

        array<int, 2> octdeg = this->mtof_01_innerScala_octdegree(degree, this->mtof_01_innerScala_kbmSize);
        int octave = (int)(octdeg[0]);
        Index index = (Index)(octdeg[1]);

        if (this->mtof_01_innerScala_kbmMapSize > index) {
            degree = (int)(this->mtof_01_innerScala_kbmValid[(Index)(this->mtof_01_innerScala_KBM_MAP_OFFSET + index)]);

            if (degree >= 0) {
                return {degree + octave * this->mtof_01_innerScala_kbmOctaveDegree, 1};
            }
        }
    }

    return {-1, 0};
}

number mtof_01_innerScala_applySCL(int degree, number frac, number refFreq) {
    array<int, 2> octdeg = this->mtof_01_innerScala_octdegree(degree, this->mtof_01_innerScala_sclEntryCount);
    return this->mtof_01_innerScala_applySCLOctIndex(octdeg[0], octdeg[1], frac, refFreq);
}

number mtof_01_innerScala_applySCLOctIndex(number octave, int index, number frac, number refFreq) {
    number p = 0;

    if (index > 0) {
        p = this->mtof_01_innerScala_sclExpMul[(Index)(index - 1)];
    }

    if (frac > 0) {
        p = this->linearinterp(frac, p, this->mtof_01_innerScala_sclExpMul[(Index)index]);
    } else if (frac < 0) {
        p = this->linearinterp(-frac, this->mtof_01_innerScala_sclExpMul[(Index)index], p);
    }

    return refFreq * rnbo_pow(2, p + octave * this->mtof_01_innerScala_sclOctaveMul);
}

void mtof_01_innerScala_updateRefFreq() {
    this->mtof_01_innerScala_lastValid = false;
    int refOffset = (int)(this->mtof_01_innerScala_kbmRefNum - this->mtof_01_innerScala_kbmMid);

    if (refOffset == 0) {
        this->mtof_01_innerScala_refFreq = this->mtof_01_innerScala_kbmRefFreq;
    } else {
        int base = (int)(this->mtof_01_innerScala_kbmSize);

        if (base < 1) {
            base = this->mtof_01_innerScala_sclEntryCount;
        }

        array<int, 2> octdeg = this->mtof_01_innerScala_octdegree(refOffset, base);
        number oct = (number)(octdeg[0]);
        int index = (int)(octdeg[1]);

        if (base > 0) {
            oct = oct + rnbo_floor(index / base);
            index = index % base;
        }

        if (index >= 0 && index < this->mtof_01_innerScala_kbmSize) {
            if (index < this->mtof_01_innerScala_kbmMapSize) {
                index = (int)(this->mtof_01_innerScala_kbmValid[(Index)((Index)(index) + this->mtof_01_innerScala_KBM_MAP_OFFSET)]);
            } else {
                index = -1;
            }
        }

        if (index < 0 || index > this->mtof_01_innerScala_sclExpMul->length)
            {} else {
            number p = 0;

            if (index > 0) {
                p = this->mtof_01_innerScala_sclExpMul[(Index)(index - 1)];
            }

            this->mtof_01_innerScala_refFreq = this->mtof_01_innerScala_kbmRefFreq / rnbo_pow(2, p + oct * this->mtof_01_innerScala_sclOctaveMul);
        }
    }
}

void mtof_01_innerScala_reset() {
    this->mtof_01_innerScala_internal = true;
    this->mtof_01_innerScala_lastValid = false;
    this->mtof_01_innerScala_lastNote = 0;
    this->mtof_01_innerScala_lastFreq = 0;
    this->mtof_01_innerScala_sclEntryCount = 0;
    this->mtof_01_innerScala_sclOctaveMul = 1;
    this->mtof_01_innerScala_sclExpMul = {};
    this->mtof_01_innerScala_kbmValid = {0, 0, 0, 60, 69, 440};
    this->mtof_01_innerScala_kbmMid = 60;
    this->mtof_01_innerScala_kbmRefNum = 69;
    this->mtof_01_innerScala_kbmRefFreq = 440;
    this->mtof_01_innerScala_kbmSize = 0;
    this->mtof_01_innerScala_kbmMin = 0;
    this->mtof_01_innerScala_kbmMax = 0;
    this->mtof_01_innerScala_kbmOctaveDegree = 12;
    this->mtof_01_innerScala_kbmMapSize = 0;
    this->mtof_01_innerScala_refFreq = 261.63;
}

void mtof_01_init() {
    this->mtof_01_innerScala_update(this->mtof_01_scale, this->mtof_01_map);
}

number ftom_01_innerFtoM_next(number frequency, number tuning) {
    if (frequency <= 0.0) {
        return -999;
    }

    if (frequency == this->ftom_01_innerFtoM_lastInValue && tuning == this->ftom_01_innerFtoM_lastTuning) {
        return this->ftom_01_innerFtoM_lastOutValue;
    }

    this->ftom_01_innerFtoM_lastInValue = frequency;
    this->ftom_01_innerFtoM_lastTuning = tuning;
    this->ftom_01_innerFtoM_lastOutValue = (frequency == 0 || tuning == 0 ? 0 : 69. + 17.31234050465299 * rnbo_log(frequency / tuning));

    {
        this->ftom_01_innerFtoM_lastOutValue = rnbo_fround(this->ftom_01_innerFtoM_lastOutValue * 1 / (number)1) * 1;
    }

    return this->ftom_01_innerFtoM_lastOutValue;
}

void ftom_01_innerFtoM_reset() {
    this->ftom_01_innerFtoM_lastInValue = 0;
    this->ftom_01_innerFtoM_lastOutValue = 0;
    this->ftom_01_innerFtoM_lastTuning = 0;
}

void ftom_01_innerScala_mid(int v) {
    this->ftom_01_innerScala_kbmMid = v;
    this->ftom_01_innerScala_updateRefFreq();
}

void ftom_01_innerScala_ref(int v) {
    this->ftom_01_innerScala_kbmRefNum = v;
    this->ftom_01_innerScala_updateRefFreq();
}

void ftom_01_innerScala_base(number v) {
    this->ftom_01_innerScala_kbmRefFreq = v;
    this->ftom_01_innerScala_updateRefFreq();
}

void ftom_01_innerScala_init() {
    list sclValid = {
        12,
        100,
        0,
        200,
        0,
        300,
        0,
        400,
        0,
        500,
        0,
        600,
        0,
        700,
        0,
        800,
        0,
        900,
        0,
        1000,
        0,
        1100,
        0,
        2,
        1
    };

    this->ftom_01_innerScala_updateScale(sclValid);
}

void ftom_01_innerScala_update(list scale, list map) {
    if (scale->length > 0) {
        this->ftom_01_innerScala_updateScale(scale);
    }

    if (map->length > 0) {
        this->ftom_01_innerScala_updateMap(map);
    }
}

number ftom_01_innerScala_mtof(number note) {
    if ((bool)(this->ftom_01_innerScala_lastValid) && this->ftom_01_innerScala_lastNote == note) {
        return this->ftom_01_innerScala_lastFreq;
    }

    array<int, 2> degoct = this->ftom_01_innerScala_applyKBM(note);
    number out = 0;

    if (degoct[1] > 0) {
        out = this->ftom_01_innerScala_applySCL(degoct[0], fract(note), this->ftom_01_innerScala_refFreq);
    }

    this->ftom_01_innerScala_updateLast(note, out);
    return out;
}

number ftom_01_innerScala_ftom(number hz) {
    if (hz <= 0.0) {
        return 0.0;
    }

    if ((bool)(this->ftom_01_innerScala_lastValid) && this->ftom_01_innerScala_lastFreq == hz) {
        return this->ftom_01_innerScala_lastNote;
    }

    array<number, 2> df = this->ftom_01_innerScala_hztodeg(hz);
    int degree = (int)(df[0]);
    number frac = df[1];
    number out = 0;

    if (this->ftom_01_innerScala_kbmSize == 0) {
        out = this->ftom_01_innerScala_kbmMid + degree;
    } else {
        array<int, 2> octdeg = this->ftom_01_innerScala_octdegree(degree, this->ftom_01_innerScala_kbmOctaveDegree);
        number oct = (number)(octdeg[0]);
        int index = (int)(octdeg[1]);
        Index entry = 0;

        for (Index i = 0; i < this->ftom_01_innerScala_kbmMapSize; i++) {
            if (index == this->ftom_01_innerScala_kbmValid[(Index)(i + this->ftom_01_innerScala_KBM_MAP_OFFSET)]) {
                entry = i;
                break;
            }
        }

        out = oct * this->ftom_01_innerScala_kbmSize + entry + this->ftom_01_innerScala_kbmMid;
    }

    out = out + frac;
    this->ftom_01_innerScala_updateLast(out, hz);
    return this->ftom_01_innerScala_lastNote;
}

int ftom_01_innerScala_updateScale(list scl) {
    if (scl->length > 1 && scl[0] * 2 + 1 == scl->length) {
        this->ftom_01_innerScala_lastValid = false;
        this->ftom_01_innerScala_sclExpMul = {};
        number last = 1;

        for (Index i = 1; i < scl->length; i += 2) {
            const number c = (const number)(scl[(Index)(i + 0)]);
            const number d = (const number)(scl[(Index)(i + 1)]);

            if (d <= 0) {
                last = c / (number)1200;
            } else {
                last = rnbo_log2(c / d);
            }

            this->ftom_01_innerScala_sclExpMul->push(last);
        }

        this->ftom_01_innerScala_sclOctaveMul = last;
        this->ftom_01_innerScala_sclEntryCount = (int)(this->ftom_01_innerScala_sclExpMul->length);
        this->ftom_01_innerScala_updateRefFreq();
        return 1;
    }

    return 0;
}

int ftom_01_innerScala_updateMap(list kbm) {
    if (kbm->length == 1 && kbm[0] == 0.0) {
        kbm = {0.0, 0.0, 0.0, 60.0, 69.0, 440.0};
    }

    if (kbm->length >= 6 && kbm[0] >= 0.0) {
        this->ftom_01_innerScala_lastValid = false;
        Index size = (Index)(kbm[0]);
        int octave = 12;

        if (kbm->length > 6) {
            octave = (int)(kbm[6]);
        }

        if (size > 0 && kbm->length < this->ftom_01_innerScala_KBM_MAP_OFFSET) {
            return 0;
        }

        this->ftom_01_innerScala_kbmSize = (int)(size);
        this->ftom_01_innerScala_kbmMin = (int)(kbm[1]);
        this->ftom_01_innerScala_kbmMax = (int)(kbm[2]);
        this->ftom_01_innerScala_kbmMid = (int)(kbm[3]);
        this->ftom_01_innerScala_kbmRefNum = (int)(kbm[4]);
        this->ftom_01_innerScala_kbmRefFreq = kbm[5];
        this->ftom_01_innerScala_kbmOctaveDegree = octave;
        this->ftom_01_innerScala_kbmValid = kbm;
        this->ftom_01_innerScala_kbmMapSize = (kbm->length - this->ftom_01_innerScala_KBM_MAP_OFFSET > kbm->length ? kbm->length : (kbm->length - this->ftom_01_innerScala_KBM_MAP_OFFSET < 0 ? 0 : kbm->length - this->ftom_01_innerScala_KBM_MAP_OFFSET));
        this->ftom_01_innerScala_updateRefFreq();
        return 1;
    }

    return 0;
}

void ftom_01_innerScala_updateLast(number note, number freq) {
    this->ftom_01_innerScala_lastValid = true;
    this->ftom_01_innerScala_lastNote = note;
    this->ftom_01_innerScala_lastFreq = freq;
}

array<number, 2> ftom_01_innerScala_hztodeg(number hz) {
    number hza = rnbo_abs(hz);

    number octave = rnbo_floor(
        rnbo_log2(hza / this->ftom_01_innerScala_refFreq) / this->ftom_01_innerScala_sclOctaveMul
    );

    int i = 0;
    number frac = 0;
    number n = 0;

    for (; i < this->ftom_01_innerScala_sclEntryCount; i++) {
        number c = this->ftom_01_innerScala_applySCLOctIndex(octave, i + 0, 0.0, this->ftom_01_innerScala_refFreq);
        n = this->ftom_01_innerScala_applySCLOctIndex(octave, i + 1, 0.0, this->ftom_01_innerScala_refFreq);

        if (c <= hza && hza < n) {
            if (c != hza) {
                frac = rnbo_log2(hza / c) / rnbo_log2(n / c);
            }

            break;
        }
    }

    if (i == this->ftom_01_innerScala_sclEntryCount && n != hza) {
        number c = n;
        n = this->ftom_01_innerScala_applySCLOctIndex(octave + 1, 0, 0.0, this->ftom_01_innerScala_refFreq);
        frac = rnbo_log2(hza / c) / rnbo_log2(n / c);
    }

    number deg = i + octave * this->ftom_01_innerScala_sclEntryCount;

    {
        deg = rnbo_fround((deg + frac) * 1 / (number)1) * 1;
        frac = 0.0;
    }

    return {deg, frac};
}

array<int, 2> ftom_01_innerScala_octdegree(int degree, int count) {
    int octave = 0;
    int index = 0;

    if (degree < 0) {
        octave = -(1 + (-1 - degree) / count);
        index = -degree % count;

        if (index > 0) {
            index = count - index;
        }
    } else {
        octave = degree / count;
        index = degree % count;
    }

    return {octave, index};
}

array<int, 2> ftom_01_innerScala_applyKBM(number note) {
    if ((this->ftom_01_innerScala_kbmMin == this->ftom_01_innerScala_kbmMax && this->ftom_01_innerScala_kbmMax == 0) || (note >= this->ftom_01_innerScala_kbmMin && note <= this->ftom_01_innerScala_kbmMax)) {
        int degree = (int)(rnbo_floor(note - this->ftom_01_innerScala_kbmMid));

        if (this->ftom_01_innerScala_kbmSize == 0) {
            return {degree, 1};
        }

        array<int, 2> octdeg = this->ftom_01_innerScala_octdegree(degree, this->ftom_01_innerScala_kbmSize);
        int octave = (int)(octdeg[0]);
        Index index = (Index)(octdeg[1]);

        if (this->ftom_01_innerScala_kbmMapSize > index) {
            degree = (int)(this->ftom_01_innerScala_kbmValid[(Index)(this->ftom_01_innerScala_KBM_MAP_OFFSET + index)]);

            if (degree >= 0) {
                return {degree + octave * this->ftom_01_innerScala_kbmOctaveDegree, 1};
            }
        }
    }

    return {-1, 0};
}

number ftom_01_innerScala_applySCL(int degree, number frac, number refFreq) {
    array<int, 2> octdeg = this->ftom_01_innerScala_octdegree(degree, this->ftom_01_innerScala_sclEntryCount);
    return this->ftom_01_innerScala_applySCLOctIndex(octdeg[0], octdeg[1], frac, refFreq);
}

number ftom_01_innerScala_applySCLOctIndex(number octave, int index, number frac, number refFreq) {
    number p = 0;

    if (index > 0) {
        p = this->ftom_01_innerScala_sclExpMul[(Index)(index - 1)];
    }

    if (frac > 0) {
        p = this->linearinterp(frac, p, this->ftom_01_innerScala_sclExpMul[(Index)index]);
    } else if (frac < 0) {
        p = this->linearinterp(-frac, this->ftom_01_innerScala_sclExpMul[(Index)index], p);
    }

    return refFreq * rnbo_pow(2, p + octave * this->ftom_01_innerScala_sclOctaveMul);
}

void ftom_01_innerScala_updateRefFreq() {
    this->ftom_01_innerScala_lastValid = false;
    int refOffset = (int)(this->ftom_01_innerScala_kbmRefNum - this->ftom_01_innerScala_kbmMid);

    if (refOffset == 0) {
        this->ftom_01_innerScala_refFreq = this->ftom_01_innerScala_kbmRefFreq;
    } else {
        int base = (int)(this->ftom_01_innerScala_kbmSize);

        if (base < 1) {
            base = this->ftom_01_innerScala_sclEntryCount;
        }

        array<int, 2> octdeg = this->ftom_01_innerScala_octdegree(refOffset, base);
        number oct = (number)(octdeg[0]);
        int index = (int)(octdeg[1]);

        if (base > 0) {
            oct = oct + rnbo_floor(index / base);
            index = index % base;
        }

        if (index >= 0 && index < this->ftom_01_innerScala_kbmSize) {
            if (index < this->ftom_01_innerScala_kbmMapSize) {
                index = (int)(this->ftom_01_innerScala_kbmValid[(Index)((Index)(index) + this->ftom_01_innerScala_KBM_MAP_OFFSET)]);
            } else {
                index = -1;
            }
        }

        if (index < 0 || index > this->ftom_01_innerScala_sclExpMul->length)
            {} else {
            number p = 0;

            if (index > 0) {
                p = this->ftom_01_innerScala_sclExpMul[(Index)(index - 1)];
            }

            this->ftom_01_innerScala_refFreq = this->ftom_01_innerScala_kbmRefFreq / rnbo_pow(2, p + oct * this->ftom_01_innerScala_sclOctaveMul);
        }
    }
}

void ftom_01_innerScala_reset() {
    this->ftom_01_innerScala_internal = true;
    this->ftom_01_innerScala_lastValid = false;
    this->ftom_01_innerScala_lastNote = 0;
    this->ftom_01_innerScala_lastFreq = 0;
    this->ftom_01_innerScala_sclEntryCount = 0;
    this->ftom_01_innerScala_sclOctaveMul = 1;
    this->ftom_01_innerScala_sclExpMul = {};
    this->ftom_01_innerScala_kbmValid = {0, 0, 0, 60, 69, 440};
    this->ftom_01_innerScala_kbmMid = 60;
    this->ftom_01_innerScala_kbmRefNum = 69;
    this->ftom_01_innerScala_kbmRefFreq = 440;
    this->ftom_01_innerScala_kbmSize = 0;
    this->ftom_01_innerScala_kbmMin = 0;
    this->ftom_01_innerScala_kbmMax = 0;
    this->ftom_01_innerScala_kbmOctaveDegree = 12;
    this->ftom_01_innerScala_kbmMapSize = 0;
    this->ftom_01_innerScala_refFreq = 261.63;
}

void ftom_01_init() {
    this->ftom_01_innerScala_update(this->ftom_01_scale, this->ftom_01_map);
}

number cycle_tilde_01_ph_next(number freq, number reset) {
    {
        {
            if (reset >= 0.)
                this->cycle_tilde_01_ph_currentPhase = reset;
        }
    }

    number pincr = freq * this->cycle_tilde_01_ph_conv;

    if (this->cycle_tilde_01_ph_currentPhase < 0.)
        this->cycle_tilde_01_ph_currentPhase = 1. + this->cycle_tilde_01_ph_currentPhase;

    if (this->cycle_tilde_01_ph_currentPhase > 1.)
        this->cycle_tilde_01_ph_currentPhase = this->cycle_tilde_01_ph_currentPhase - 1.;

    number tmp = this->cycle_tilde_01_ph_currentPhase;
    this->cycle_tilde_01_ph_currentPhase += pincr;
    return tmp;
}

void cycle_tilde_01_ph_reset() {
    this->cycle_tilde_01_ph_currentPhase = 0;
}

void cycle_tilde_01_ph_dspsetup() {
    this->cycle_tilde_01_ph_conv = (number)1 / this->sr;
}

void cycle_tilde_01_dspsetup(bool force) {
    if ((bool)(this->cycle_tilde_01_setupDone) && (bool)(!(bool)(force)))
        return;

    this->cycle_tilde_01_phasei = 0;
    this->cycle_tilde_01_f2i = (number)4294967296 / this->samplerate();
    this->cycle_tilde_01_wrap = (long)(this->cycle_tilde_01_buffer->getSize()) - 1;
    this->cycle_tilde_01_setupDone = true;
    this->cycle_tilde_01_ph_dspsetup();
}

void cycle_tilde_01_bufferUpdated() {
    this->cycle_tilde_01_wrap = (long)(this->cycle_tilde_01_buffer->getSize()) - 1;
}

void param_01_getPresetValue(PatcherStateInterface& preset) {
    preset["value"] = this->param_01_value;
}

void param_01_setPresetValue(PatcherStateInterface& preset) {
    if ((bool)(stateIsEmpty(preset)))
        return;

    this->param_01_value_set(preset["value"]);
}

number globaltransport_getTempoAtSample(SampleIndex sampleOffset) {
    RNBO_UNUSED(sampleOffset);
    return (this->vs > 0 ? this->globaltransport_tempo[(Index)sampleOffset] : this->globaltransport_lastTempo);
}

number globaltransport_getTempo() {
    return this->globaltransport_getTempoAtSample(this->sampleOffsetIntoNextAudioBuffer);
}

number globaltransport_getStateAtSample(SampleIndex sampleOffset) {
    RNBO_UNUSED(sampleOffset);
    return (this->vs > 0 ? this->globaltransport_state[(Index)sampleOffset] : this->globaltransport_lastState);
}

number globaltransport_getState() {
    return this->globaltransport_getStateAtSample(this->sampleOffsetIntoNextAudioBuffer);
}

number globaltransport_getBeatTimeAtMsTime(MillisecondTime time) {
    number i = 2;

    while (i < this->globaltransport_beatTimeChanges->length && this->globaltransport_beatTimeChanges[(Index)(i + 1)] <= time) {
        i += 2;
    }

    i -= 2;
    number beatTimeBase = this->globaltransport_beatTimeChanges[(Index)i];

    if (this->globaltransport_getState() == 0)
        return beatTimeBase;

    number beatTimeBaseMsTime = this->globaltransport_beatTimeChanges[(Index)(i + 1)];
    number diff = time - beatTimeBaseMsTime;
    return beatTimeBase + this->mstobeats(diff);
}

bool globaltransport_setTempo(number tempo, bool notify) {
    if ((bool)(notify)) {
        this->processTempoEvent(this->currenttime(), tempo);
        this->globaltransport_notify = true;
    } else if (this->globaltransport_getTempo() != tempo) {
        MillisecondTime ct = this->currenttime();
        this->globaltransport_beatTimeChanges->push(this->globaltransport_getBeatTimeAtMsTime(ct));
        this->globaltransport_beatTimeChanges->push(ct);

        fillSignal(
            this->globaltransport_tempo,
            this->vs,
            tempo,
            (Index)(this->sampleOffsetIntoNextAudioBuffer)
        );

        this->globaltransport_lastTempo = tempo;
        this->globaltransport_tempoNeedsReset = true;
        return true;
    }

    return false;
}

number globaltransport_getBeatTime() {
    return this->globaltransport_getBeatTimeAtMsTime(this->currenttime());
}

bool globaltransport_setState(number state, bool notify) {
    if ((bool)(notify)) {
        this->processTransportEvent(this->currenttime(), TransportState(state));
        this->globaltransport_notify = true;
    } else if (this->globaltransport_getState() != state) {
        fillSignal(
            this->globaltransport_state,
            this->vs,
            state,
            (Index)(this->sampleOffsetIntoNextAudioBuffer)
        );

        this->globaltransport_lastState = TransportState(state);
        this->globaltransport_stateNeedsReset = true;

        if (state == 0) {
            this->globaltransport_beatTimeChanges->push(this->globaltransport_getBeatTime());
            this->globaltransport_beatTimeChanges->push(this->currenttime());
        }

        return true;
    }

    return false;
}

bool globaltransport_setBeatTime(number beattime, bool notify) {
    if ((bool)(notify)) {
        this->processBeatTimeEvent(this->currenttime(), beattime);
        this->globaltransport_notify = true;
        return false;
    } else {
        bool beatTimeHasChanged = false;
        float oldBeatTime = (float)(this->globaltransport_getBeatTime());
        float newBeatTime = (float)(beattime);

        if (oldBeatTime != newBeatTime) {
            beatTimeHasChanged = true;
        }

        this->globaltransport_beatTimeChanges->push(beattime);
        this->globaltransport_beatTimeChanges->push(this->currenttime());
        return beatTimeHasChanged;
    }
}

number globaltransport_getBeatTimeAtSample(SampleIndex sampleOffset) {
    MillisecondTime msOffset = this->sampstoms(sampleOffset);
    return this->globaltransport_getBeatTimeAtMsTime(this->currenttime() + msOffset);
}

array<number, 2> globaltransport_getTimeSignatureAtMsTime(MillisecondTime time) {
    number i = 3;

    while (i < this->globaltransport_timeSignatureChanges->length && this->globaltransport_timeSignatureChanges[(Index)(i + 2)] <= time) {
        i += 3;
    }

    i -= 3;

    return {
        this->globaltransport_timeSignatureChanges[(Index)i],
        this->globaltransport_timeSignatureChanges[(Index)(i + 1)]
    };
}

array<number, 2> globaltransport_getTimeSignature() {
    return this->globaltransport_getTimeSignatureAtMsTime(this->currenttime());
}

array<number, 2> globaltransport_getTimeSignatureAtSample(SampleIndex sampleOffset) {
    MillisecondTime msOffset = this->sampstoms(sampleOffset);
    return this->globaltransport_getTimeSignatureAtMsTime(this->currenttime() + msOffset);
}

bool globaltransport_setTimeSignature(number numerator, number denominator, bool notify) {
    if ((bool)(notify)) {
        this->processTimeSignatureEvent(this->currenttime(), (int)(numerator), (int)(denominator));
        this->globaltransport_notify = true;
    } else {
        array<number, 2> currentSig = this->globaltransport_getTimeSignature();

        if (currentSig[0] != numerator || currentSig[1] != denominator) {
            this->globaltransport_timeSignatureChanges->push(numerator);
            this->globaltransport_timeSignatureChanges->push(denominator);
            this->globaltransport_timeSignatureChanges->push(this->currenttime());
            return true;
        }
    }

    return false;
}

void globaltransport_advance() {
    if ((bool)(this->globaltransport_tempoNeedsReset)) {
        fillSignal(this->globaltransport_tempo, this->vs, this->globaltransport_lastTempo);
        this->globaltransport_tempoNeedsReset = false;

        if ((bool)(this->globaltransport_notify)) {
            this->getEngine()->sendTempoEvent(this->globaltransport_lastTempo);
        }
    }

    if ((bool)(this->globaltransport_stateNeedsReset)) {
        fillSignal(this->globaltransport_state, this->vs, this->globaltransport_lastState);
        this->globaltransport_stateNeedsReset = false;

        if ((bool)(this->globaltransport_notify)) {
            this->getEngine()->sendTransportEvent(TransportState(this->globaltransport_lastState));
        }
    }

    if (this->globaltransport_beatTimeChanges->length > 2) {
        this->globaltransport_beatTimeChanges[0] = this->globaltransport_beatTimeChanges[(Index)(this->globaltransport_beatTimeChanges->length - 2)];
        this->globaltransport_beatTimeChanges[1] = this->globaltransport_beatTimeChanges[(Index)(this->globaltransport_beatTimeChanges->length - 1)];
        this->globaltransport_beatTimeChanges->length = 2;

        if ((bool)(this->globaltransport_notify)) {
            this->getEngine()->sendBeatTimeEvent(this->globaltransport_beatTimeChanges[0]);
        }
    }

    if (this->globaltransport_timeSignatureChanges->length > 3) {
        this->globaltransport_timeSignatureChanges[0] = this->globaltransport_timeSignatureChanges[(Index)(this->globaltransport_timeSignatureChanges->length - 3)];
        this->globaltransport_timeSignatureChanges[1] = this->globaltransport_timeSignatureChanges[(Index)(this->globaltransport_timeSignatureChanges->length - 2)];
        this->globaltransport_timeSignatureChanges[2] = this->globaltransport_timeSignatureChanges[(Index)(this->globaltransport_timeSignatureChanges->length - 1)];
        this->globaltransport_timeSignatureChanges->length = 3;

        if ((bool)(this->globaltransport_notify)) {
            this->getEngine()->sendTimeSignatureEvent(
                (int)(this->globaltransport_timeSignatureChanges[0]),
                (int)(this->globaltransport_timeSignatureChanges[1])
            );
        }
    }

    this->globaltransport_notify = false;
}

void globaltransport_dspsetup(bool force) {
    if ((bool)(this->globaltransport_setupDone) && (bool)(!(bool)(force)))
        return;

    fillSignal(this->globaltransport_tempo, this->vs, this->globaltransport_lastTempo);
    this->globaltransport_tempoNeedsReset = false;
    fillSignal(this->globaltransport_state, this->vs, this->globaltransport_lastState);
    this->globaltransport_stateNeedsReset = false;
    this->globaltransport_setupDone = true;
}

bool stackprotect_check() {
    this->stackprotect_count++;

    if (this->stackprotect_count > 128) {
        console->log("STACK OVERFLOW DETECTED - stopped processing branch !");
        return true;
    }

    return false;
}

void updateTime(MillisecondTime time) {
    this->_currentTime = time;
    this->sampleOffsetIntoNextAudioBuffer = (SampleIndex)(rnbo_fround(this->msToSamps(time - this->getEngine()->getCurrentTime(), this->sr)));

    if (this->sampleOffsetIntoNextAudioBuffer >= (SampleIndex)(this->vs))
        this->sampleOffsetIntoNextAudioBuffer = (SampleIndex)(this->vs) - 1;

    if (this->sampleOffsetIntoNextAudioBuffer < 0)
        this->sampleOffsetIntoNextAudioBuffer = 0;
}

void assign_defaults()
{
    mtof_01_base = 440;
    ftom_01_base = 440;
    kink_tilde_01_x = 0;
    kink_tilde_01_slope = 0.5;
    cycle_tilde_01_frequency = 0;
    cycle_tilde_01_phase_offset = 0;
    param_01_value = 220;
    _currentTime = 0;
    audioProcessSampleCount = 0;
    sampleOffsetIntoNextAudioBuffer = 0;
    zeroBuffer = nullptr;
    dummyBuffer = nullptr;
    signals[0] = nullptr;
    didAllocateSignals = 0;
    vs = 0;
    maxvs = 0;
    sr = 44100;
    invsr = 0.00002267573696;
    mtof_01_innerMtoF_lastInValue = 0;
    mtof_01_innerMtoF_lastOutValue = 0;
    mtof_01_innerMtoF_lastTuning = 0;
    mtof_01_innerScala_internal = true;
    mtof_01_innerScala_lastValid = false;
    mtof_01_innerScala_lastNote = 0;
    mtof_01_innerScala_lastFreq = 0;
    mtof_01_innerScala_sclEntryCount = 0;
    mtof_01_innerScala_sclOctaveMul = 1;
    mtof_01_innerScala_kbmValid = { 0, 0, 0, 60, 69, 440 };
    mtof_01_innerScala_kbmMid = 60;
    mtof_01_innerScala_kbmRefNum = 69;
    mtof_01_innerScala_kbmRefFreq = 440;
    mtof_01_innerScala_kbmSize = 0;
    mtof_01_innerScala_kbmMin = 0;
    mtof_01_innerScala_kbmMax = 0;
    mtof_01_innerScala_kbmOctaveDegree = 12;
    mtof_01_innerScala_kbmMapSize = 0;
    mtof_01_innerScala_refFreq = 261.63;
    ftom_01_innerFtoM_lastInValue = 0;
    ftom_01_innerFtoM_lastOutValue = 0;
    ftom_01_innerFtoM_lastTuning = 0;
    ftom_01_innerScala_internal = true;
    ftom_01_innerScala_lastValid = false;
    ftom_01_innerScala_lastNote = 0;
    ftom_01_innerScala_lastFreq = 0;
    ftom_01_innerScala_sclEntryCount = 0;
    ftom_01_innerScala_sclOctaveMul = 1;
    ftom_01_innerScala_kbmValid = { 0, 0, 0, 60, 69, 440 };
    ftom_01_innerScala_kbmMid = 60;
    ftom_01_innerScala_kbmRefNum = 69;
    ftom_01_innerScala_kbmRefFreq = 440;
    ftom_01_innerScala_kbmSize = 0;
    ftom_01_innerScala_kbmMin = 0;
    ftom_01_innerScala_kbmMax = 0;
    ftom_01_innerScala_kbmOctaveDegree = 12;
    ftom_01_innerScala_kbmMapSize = 0;
    ftom_01_innerScala_refFreq = 261.63;
    cycle_tilde_01_wrap = 0;
    cycle_tilde_01_ph_currentPhase = 0;
    cycle_tilde_01_ph_conv = 0;
    cycle_tilde_01_setupDone = false;
    param_01_lastValue = 0;
    globaltransport_tempo = nullptr;
    globaltransport_tempoNeedsReset = false;
    globaltransport_lastTempo = 120;
    globaltransport_state = nullptr;
    globaltransport_stateNeedsReset = false;
    globaltransport_lastState = 0;
    globaltransport_beatTimeChanges = { 0, 0 };
    globaltransport_timeSignatureChanges = { 4, 4, 0 };
    globaltransport_notify = false;
    globaltransport_setupDone = false;
    stackprotect_count = 0;
    _voiceIndex = 0;
    _noteNumber = 0;
    isMuted = 1;
}

// member variables

    list mtof_01_midivalue;
    list mtof_01_scale;
    list mtof_01_map;
    number mtof_01_base;
    list ftom_01_frequency;
    list ftom_01_scale;
    list ftom_01_map;
    number ftom_01_base;
    number kink_tilde_01_x;
    number kink_tilde_01_slope;
    number cycle_tilde_01_frequency;
    number cycle_tilde_01_phase_offset;
    number param_01_value;
    MillisecondTime _currentTime;
    SampleIndex audioProcessSampleCount;
    SampleIndex sampleOffsetIntoNextAudioBuffer;
    signal zeroBuffer;
    signal dummyBuffer;
    SampleValue * signals[1];
    bool didAllocateSignals;
    Index vs;
    Index maxvs;
    number sr;
    number invsr;
    number mtof_01_innerMtoF_lastInValue;
    number mtof_01_innerMtoF_lastOutValue;
    number mtof_01_innerMtoF_lastTuning;
    Float64BufferRef mtof_01_innerMtoF_buffer;
    bool mtof_01_innerScala_internal;
    const Index mtof_01_innerScala_KBM_MAP_OFFSET = 7;
    bool mtof_01_innerScala_lastValid;
    number mtof_01_innerScala_lastNote;
    number mtof_01_innerScala_lastFreq;
    int mtof_01_innerScala_sclEntryCount;
    number mtof_01_innerScala_sclOctaveMul;
    list mtof_01_innerScala_sclExpMul;
    list mtof_01_innerScala_kbmValid;
    int mtof_01_innerScala_kbmMid;
    int mtof_01_innerScala_kbmRefNum;
    number mtof_01_innerScala_kbmRefFreq;
    int mtof_01_innerScala_kbmSize;
    int mtof_01_innerScala_kbmMin;
    int mtof_01_innerScala_kbmMax;
    int mtof_01_innerScala_kbmOctaveDegree;
    Index mtof_01_innerScala_kbmMapSize;
    number mtof_01_innerScala_refFreq;
    number ftom_01_innerFtoM_lastInValue;
    number ftom_01_innerFtoM_lastOutValue;
    number ftom_01_innerFtoM_lastTuning;
    bool ftom_01_innerScala_internal;
    const Index ftom_01_innerScala_KBM_MAP_OFFSET = 7;
    bool ftom_01_innerScala_lastValid;
    number ftom_01_innerScala_lastNote;
    number ftom_01_innerScala_lastFreq;
    int ftom_01_innerScala_sclEntryCount;
    number ftom_01_innerScala_sclOctaveMul;
    list ftom_01_innerScala_sclExpMul;
    list ftom_01_innerScala_kbmValid;
    int ftom_01_innerScala_kbmMid;
    int ftom_01_innerScala_kbmRefNum;
    number ftom_01_innerScala_kbmRefFreq;
    int ftom_01_innerScala_kbmSize;
    int ftom_01_innerScala_kbmMin;
    int ftom_01_innerScala_kbmMax;
    int ftom_01_innerScala_kbmOctaveDegree;
    Index ftom_01_innerScala_kbmMapSize;
    number ftom_01_innerScala_refFreq;
    Float64BufferRef cycle_tilde_01_buffer;
    long cycle_tilde_01_wrap;
    uint32_t cycle_tilde_01_phasei;
    SampleValue cycle_tilde_01_f2i;
    number cycle_tilde_01_ph_currentPhase;
    number cycle_tilde_01_ph_conv;
    bool cycle_tilde_01_setupDone;
    number param_01_lastValue;
    signal globaltransport_tempo;
    bool globaltransport_tempoNeedsReset;
    number globaltransport_lastTempo;
    signal globaltransport_state;
    bool globaltransport_stateNeedsReset;
    number globaltransport_lastState;
    list globaltransport_beatTimeChanges;
    list globaltransport_timeSignatureChanges;
    bool globaltransport_notify;
    bool globaltransport_setupDone;
    number stackprotect_count;
    DataRef RNBODefaultMtofLookupTable256;
    DataRef RNBODefaultSinus;
    Index _voiceIndex;
    Int _noteNumber;
    Index isMuted;
    indexlist paramInitIndices;
    indexlist paramInitOrder;

};

PatcherInterface* creaternbomatic()
{
    return new rnbomatic();
}

#ifndef RNBO_NO_PATCHERFACTORY

extern "C" PatcherFactoryFunctionPtr GetPatcherFactoryFunction(PlatformInterface* platformInterface)
#else

extern "C" PatcherFactoryFunctionPtr rnbomaticFactoryFunction(PlatformInterface* platformInterface)
#endif

{
    Platform::set(platformInterface);
    return creaternbomatic;
}

} // end RNBO namespace

