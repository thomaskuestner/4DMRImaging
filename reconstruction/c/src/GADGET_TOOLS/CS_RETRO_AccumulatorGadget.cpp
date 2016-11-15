/*	
file name	: 	CS_Retro.cpp
author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
version		: 	v1.0
date		: 	
description	: 	
references	:	 ...cpp (origial Gadgetron - mri_core project)
changes		:	
*/

#include "CS_Retro_AccumulatorGadget.h"
#include "GadgetIsmrmrdReadWrite.h"
#include <cmath>
namespace Gadgetron{
// class constructor
CS_Retro_AccumulatorGadget::CS_Retro_AccumulatorGadget() : bufferkSpace_(0), bufferNav_(0), iBaseRes_(0.0), fFullySa_(.065), fTR_(0.0), iEchoLine_(0.0), iEchoPartition_(0.0), iNavPeriod_(0.0), iNavPERes_(0.0), lNoScans_(0), iNoSamples_(0), iNoChannels_(0) {
	CS_GlobalVar::instance()->vPE_.clear();
	CS_GlobalVar::instance()->vPA_.clear();
}
 
// class destructor - delete temporal buffer/memory
CS_Retro_AccumulatorGadget::~CS_Retro_AccumulatorGadget(){
	if (bufferkSpace_) delete bufferkSpace_;
	if (bufferNav_) delete bufferNav_;
	//GlobalVar::instance()->vPE_.clear();
	//GlobalVar::instance()->vPA_.clear();
}

// read flexible data header
int CS_Retro_AccumulatorGadget::process_config(ACE_Message_Block* mb)
{
	// read xml header file
	boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));
	
	// create sequence encoding parameters
	ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
	ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
	ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
	ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

	// get matrix size
	dimensionsIn_.push_back(r_space.matrixSize().x());
	dimensionsIn_.push_back(e_space.matrixSize().y());
	dimensionsIn_.push_back(e_space.matrixSize().z());
	GADGET_DEBUG2("Matrix size: %d, %d, %d\n", dimensionsIn_[0], dimensionsIn_[1], dimensionsIn_[2]);
	
	// get FOV
	field_of_view_.push_back(r_space.fieldOfView_mm().x());
	field_of_view_.push_back(e_space.fieldOfView_mm().y());
	field_of_view_.push_back(e_space.fieldOfView_mm().z());
	GADGET_DEBUG2("FOV: %f, %f, %f\n", r_space.fieldOfView_mm().x(), e_space.fieldOfView_mm().y(), e_space.fieldOfView_mm().z());

	// get number of repetitions
	repetitions_ = e_limits.repetition().present() ? e_limits.repetition().get().maximum()+1 : 1;
	dimensionsIn_.push_back(repetitions_);
	GADGET_DEBUG2("Number of repetitions: %i\n", dimensionsIn_[3]);
	
	// get number of phases
	iNPhases_ = e_limits.phase().present() ? e_limits.phase().get().maximum() : 1;
	GADGET_DEBUG2("Number of phases: %i\n", iNPhases_);

	// get number of slices
	slices_ = e_limits.slice().present() ? e_limits.slice().get().maximum()+1 : 1;
	
	// get echo line and echo partition
	iEchoLine_ = e_limits.kspace_encoding_step_1().get().center();
	iEchoPartition_ = e_limits.kspace_encoding_step_2().get().center();
	GADGET_DEBUG2("echo line: %i, echo partition: %i", iEchoLine_, iEchoPartition_);

	// read Compressed Sensing and CS_Retro values
	int iESPReSSoY = 0;
	int iESPReSSoZ = 0;
	if ((*e_seq.begin()).trajectoryDescription().present()) {
		GADGET_DEBUG1("\n\nTrajectory description present!\n\n");
		ISMRMRD::trajectoryDescriptionType traj_desc = (*e_seq.begin()).trajectoryDescription().get();

		for (ISMRMRD::trajectoryDescriptionType::userParameterLong_sequence::iterator i (traj_desc.userParameterLong().begin ()); i != traj_desc.userParameterLong().end(); ++i) {
			/*if (std::strcmp(i->name().c_str(),"NPhases") == 0) {
				iNPhases_ = i->value();
				GADGET_DEBUG2("number of phases is %i \n", iNPhases_);
			} */
			//iNPhases_ = 4;

			if (std::strcmp(i->name().c_str(),"NavPeriod") == 0) {
				iNavPeriod_ = i->value();
				GADGET_DEBUG2("nav period is %i \n", iNavPeriod_);
			} 
			
			if (std::strcmp(i->name().c_str(),"NavPERes") == 0) {
				iNavPERes_ = i->value();
				GADGET_DEBUG2("nav PE res is %i \n", iNavPERes_);
			}

			if (std::strcmp(i->name().c_str(),"ESPReSSoY") == 0) {
				iESPReSSoY = i->value();
				GADGET_DEBUG2("ESPReSSoZ res is %i \n", iESPReSSoY);
			}

			if (std::strcmp(i->name().c_str(),"ESPReSSoZ") == 0) {
				iESPReSSoZ = i->value();
				GADGET_DEBUG2("ESPReSSoZ is %i \n", iESPReSSoZ);
			}

		}

		iESPReSSoDirection_ = 10;
		fPartialFourierVal_ = 1.0;
		if ((iESPReSSoY > 9 && iESPReSSoY < 14) || (iESPReSSoZ > 9 && iESPReSSoZ < 14)) {
			GADGET_DEBUG1("Partial Fourier data..\n");
			GADGET_DEBUG2("ESPReSSo Y: %f, ESPReSSo Z: %f\n", iESPReSSoY, iESPReSSoZ);
			// get Partial Fourier dimension
			if (iESPReSSoY > 9){
				iESPReSSoDirection_ = 1;
				// get Partial Fourier value
				switch (iESPReSSoY){
					case 10:
						fPartialFourierVal_ = 0.5;
						break;
					case 11:
						fPartialFourierVal_ = 0.625;
						break;
					case 12:
						fPartialFourierVal_ = 0.75;
						break;
					case 13:
						fPartialFourierVal_ = 0.875;
						break;
					default:
						fPartialFourierVal_ = 1.0;
						break;
				}
			}
			else if (iESPReSSoZ > 9){
				iESPReSSoDirection_ = 2;
				// get Partial Fourier value
				switch (iESPReSSoZ){
					case 10:
						fPartialFourierVal_ = 0.5;
						break;
					case 11:
						fPartialFourierVal_ = 0.625;
						break;
					case 12:
						fPartialFourierVal_ = 0.75;
						break;
					case 13:
						fPartialFourierVal_ = 0.875;
						break;
					default:
						fPartialFourierVal_ = 1.0;
						break;
				}
			}
		}
	}
	else{
		GADGET_DEBUG1("\n\nNo trajectory description present!\n\n");
	}

	// repetition time
	fTR_ = cfg->sequenceParameters().get().TR().at(0);

	// get population mode
	iPopulationMode_ = this->get_int_value("PopulationMode");

	// get gating mode
	iGatingMode_ = this->get_int_value("GatingMode");

	// get maximum number of stored acquisition header information in global vector
	iAcqSize_ = this->get_int_value("Acquisition Vector Size");

	return GADGET_OK;
}

// process data - incoming unordered k-space(RO,t,c) --> ordered k-space(x,y,z,Phases,c)
int CS_Retro_AccumulatorGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2){
	
	/*---------------------------------------------------*/
	/*----------- init buffer for k-space data ----------*/
	/*---------------------------------------------------*/
	if (!bufferkSpace_){

		// concat higher and lower bytes from total measurement variable
		lNoScans_ = (((m1->getObjectPtr()->user_int[3]) << 16) & 0xFFFF0000) | (m1->getObjectPtr()->user_int[2] & 0x0000FFFF);
		lNoScans_ = 2849;

		iNPhases_ = m1->getObjectPtr()->user_int[0];
		GADGET_DEBUG2("no. of gates: %i\n", iNPhases_);
		iNPhases_ = 4;
		iNoChannels_ = m1->getObjectPtr()->active_channels;

		// initialize k-space buffer
		if (!(bufferkSpace_ = new hoNDArray< std::complex<float> >())) {
			GADGET_DEBUG1("CS_Retro: Failed to create k-space buffer\n");
			return GADGET_FAIL;
		}

		// get number of samples in acquisition (equals base resolution)
		iBaseRes_ = m1->getObjectPtr()->number_of_samples;

		GADGET_DEBUG2("base res.: %f, no. scans: %i, no. channel: %i\n", iBaseRes_, lNoScans_, iNoChannels_);
		// dimension vector of k-space array
		dimkSpace_.push_back(iBaseRes_);
		dimkSpace_.push_back(lNoScans_);
		dimkSpace_.push_back(iNoChannels_);
	
		// create buffer array for incoming k-space data (readout, time, channel)
		try {bufferkSpace_->create(&dimkSpace_);}
		catch (std::runtime_error &err){
			GADGET_DEBUG_EXCEPTION(err, "CS_Retro: Failed to allocate k-space buffer array\n");
			return GADGET_FAIL;
		}

		// copy header information of first acquisition to global variable (create new header, set memory zero, copy header, push onto global vector)
		GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* tmp_m1 = new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
		memset(tmp_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));
		copy_header(tmp_m1, m1);
		CS_GlobalVar::instance()->AcqVec_.push_back(tmp_m1);
		GADGET_DEBUG1("Receiving data..\n");

		GADGET_DEBUG1("\n \n bufferkSpace_ \n \n:");
		bufferkSpace_->print(std::cout);
	}
	

	/*---------------------------------------------------*/
	/*--------- init buffer for navigator data ----------*/
	/*---------------------------------------------------*/
	if (!bufferNav_){	
		// initialize k-space buffer
		if (!(bufferNav_ = new hoNDArray< std::complex<float> >())) {
			GADGET_DEBUG1("CS_Retro: Failed to create navigator buffer\n");
			return GADGET_FAIL;
		}
		
		// dimension vector of navigator array
		dimNav_.push_back(iBaseRes_);
		dimNav_.push_back(lNoScans_);
		dimNav_.push_back(iNavPERes_);
		dimNav_.push_back(iNoChannels_);
		iNoNav_ = 0;
		iNoNavLine_ = 0;
		GADGET_DEBUG2("navigator dimensions: base res: %i, no. scans: %i, PE resolution: %i, no. channels: %i\n", iBaseRes_, lNoScans_, iNavPERes_, iNoChannels_);

		// create buffer array for incoming navigator data (readout, time, PE, channel)
		try {bufferNav_->create(&dimNav_);}
		catch (std::runtime_error &err){
			GADGET_DEBUG_EXCEPTION(err, "CS_Retro: Failed to allocate navigator buffer array\n");
			return GADGET_FAIL;
		}
		GADGET_DEBUG1("\n \n bufferNav \n \n:");
		bufferNav_->print(std::cout);
		lCurrentScan_ = 0;
	} 

	/*---------------------------------------------------*/
	/*----------- store incoming k-space data -----------*/
	/*---------------------------------------------------*/
	// navigator flag
	bool  bNavigator = false;
	if (m1->getObjectPtr()->user_int[1] & 0x1 != 0){
		bNavigator = true;
	}

	// get current loop counters
	int samples		= m1->getObjectPtr()->number_of_samples;
	int line		= m1->getObjectPtr()->idx.kspace_encode_step_1;
	int partition	= m1->getObjectPtr()->idx.kspace_encode_step_2;
	int slice		= m1->getObjectPtr()->idx.slice;		
	int repetition	= m1->getObjectPtr()->idx.repetition;
	int set			= m1->getObjectPtr()->idx.set;
	//int scan		= m1->getObjectPtr()->scan_counter;

	// push current loop counters on according vector (temporal)
	CS_GlobalVar::instance()->vPE_.push_back(line); CS_GlobalVar::instance()->vPA_.push_back(partition);

	// get data pointer
	std::complex<float> *pkSpace	= bufferkSpace_->get_data_ptr();
	std::complex<float> *pNav		= bufferNav_->get_data_ptr();
	std::complex<float>	*pIncoming	= m2->getObjectPtr()->get_data_ptr();
	for (int c = 0; c < m1->getObjectPtr()->active_channels; c++){
		size_t offset_kSpace = c*dimkSpace_[0]*dimkSpace_[1] + lCurrentScan_*dimkSpace_[0] + (dimkSpace_[0]>>1)-m1->getObjectPtr()->center_sample;
		memcpy(pkSpace + offset_kSpace, pIncoming+c*samples, sizeof(std::complex<float>)*samples);
		
		if (bNavigator == true){		
			size_t offset_Nav = c*dimNav_[0]*dimNav_[1]*dimNav_[2] + iNoNavLine_*dimNav_[0]*dimNav_[1] + iNoNav_*dimNav_[0]+(dimNav_[0]>>1)-m1->getObjectPtr()->center_sample;
			memcpy(pNav + offset_Nav, pIncoming+c*samples, sizeof(std::complex<float>)*samples);		
		}	
	}
	if (bNavigator == true){
		iNoNavLine_++;
		if (iNoNavLine_ == iNavPERes_){
				
			iNoNav_++;
			iNoNavLine_ = 0;
		}
		if (iNoNavLine_ == (int)((float)iNavPERes_/2)){
			CS_GlobalVar::instance()->vNavInd_.push_back((float)lCurrentScan_);
		}
	}
	lCurrentScan_++;

	/*---------------------------------------------------*/
	/*--------------- process sampled data --------------*/
	/*---------------------------------------------------*/
	if (lCurrentScan_ == lNoScans_){
		GADGET_DEBUG1("data received.. try to process data\n");
		// crop non-empty data from navigator array	
		//check if last measurement is in between a navigator block
		if (iNoNavLine_!=0){
			iNoNav_--;
			CS_GlobalVar::instance()->vNavInd_.pop_back();
		}
		GADGET_DEBUG2("%i navigator data found..\n", iNoNav_);
		std::vector<size_t> vStart, vSize;
		vStart.push_back(0);  vStart.push_back(0); vStart.push_back(0); vStart.push_back(0);
		vSize.push_back(bufferNav_->get_size(0)); vSize.push_back(iNoNav_); vSize.push_back(bufferNav_->get_size(2)); vSize.push_back(bufferNav_->get_size(3)); 
		bufferNav_->print(std::cout);
		get_subarray(*bufferNav_, vStart, vSize, *bufferNav_);
		bufferNav_->print(std::cout);
		//save_array(*bufferNav_, "D:/Martin/bufferNav");
		GADGET_DEBUG1("Flag - last in measurement detected..\n");
		GADGET_DEBUG1("\n------------- navigator buffer ------------------\n");
		bufferNav_->print(std::cout);
		GADGET_DEBUG1("\n------------- kSpace buffer -----------------\n");
		bufferkSpace_->print(std::cout);

		// create new ContainerMessages for header, navigator and kSpace data
		// header
		GadgetContainerMessage<ISMRMRD::ImageHeader>* tmp_m1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		
		// initialize the image header
		memset(tmp_m1->getObjectPtr(),0,sizeof(ISMRMRD::ImageHeader));

		// initialize flags
		tmp_m1->getObjectPtr()->flags = 0;

		//tmp_m1->getObjectPtr()->user_int[0] = 7;
		tmp_m1->getObjectPtr()->user_int[0]			= iNPhases_;
		tmp_m1->getObjectPtr()->user_int[1]			= m1->getObjectPtr()->user_int[1];
		tmp_m1->getObjectPtr()->user_int[2]			= m1->getObjectPtr()->user_int[2];
		tmp_m1->getObjectPtr()->user_int[3]			= m1->getObjectPtr()->user_int[3];
		tmp_m1->getObjectPtr()->matrix_size[0]     = dimensionsIn_[0];
		tmp_m1->getObjectPtr()->matrix_size[1]     = dimensionsIn_[1];
		tmp_m1->getObjectPtr()->matrix_size[2]     = dimensionsIn_[2];
		tmp_m1->getObjectPtr()->field_of_view[0]   = field_of_view_[0];
		tmp_m1->getObjectPtr()->field_of_view[1]   = field_of_view_[1];
		tmp_m1->getObjectPtr()->field_of_view[2]   = field_of_view_[2];
		tmp_m1->getObjectPtr()->channels           = (uint16_t)m1->getObjectPtr()->active_channels;
		tmp_m1->getObjectPtr()->slice				= m1->getObjectPtr()->idx.slice;
		memcpy(tmp_m1->getObjectPtr()->position,m1->getObjectPtr()->position,sizeof(float)*3);
		memcpy(tmp_m1->getObjectPtr()->read_dir,m1->getObjectPtr()->read_dir,sizeof(float)*3);
		memcpy(tmp_m1->getObjectPtr()->phase_dir,m1->getObjectPtr()->phase_dir,sizeof(float)*3);
		memcpy(tmp_m1->getObjectPtr()->slice_dir,m1->getObjectPtr()->slice_dir, sizeof(float)*3);
		memcpy(tmp_m1->getObjectPtr()->patient_table_position,m1->getObjectPtr()->patient_table_position, sizeof(float)*3);
		tmp_m1->getObjectPtr()->image_data_type = ISMRMRD::DATA_COMPLEX_FLOAT;
		tmp_m1->getObjectPtr()->image_index = (uint16_t)(++image_counter_);
		tmp_m1->getObjectPtr()->image_series_index = (uint16_t)image_series_;

		// set user values, if Compressed Sensing is active
		if(this->get_bool_value("CS_on") == true){
			tmp_m1->getObjectPtr()->user_float[0] = fCSAcc_;
			tmp_m1->getObjectPtr()->user_float[1] = fFullySa_/100;
			tmp_m1->getObjectPtr()->user_float[2] = fPartialFourierVal_;
			tmp_m1->getObjectPtr()->user_float[3] = fLQ_;
			tmp_m1->getObjectPtr()->user_float[4] = fLESPReSSo_;

			tmp_m1->getObjectPtr()->user_int[1] = iBodyRegion_;
			tmp_m1->getObjectPtr()->user_int[2] = iSamplingType_;
			tmp_m1->getObjectPtr()->user_int[3] = iVDMap_;
			tmp_m1->getObjectPtr()->user_int[4] = iESPReSSoDirection_;
		}
		tmp_m1->getObjectPtr()->user_int[5] = iNoNav_;
		// navigator
		GadgetContainerMessage<hoNDArray<std::complex<float>>>* tmp_m2 = new GadgetContainerMessage<hoNDArray<std::complex<float>>>();

		// concatenate data with header
		tmp_m1->cont(tmp_m2);
	
		// create output
		try{ tmp_m2->getObjectPtr()->create(bufferNav_->get_dimensions());}
		catch (std::runtime_error &err){
			GADGET_DEBUG_EXCEPTION(err, "CS_Retro: Unable to allocate new image array\m");
			tmp_m1->release();
			return -1;
		}

		// copy data
		memcpy(tmp_m2->getObjectPtr()->get_data_ptr(), bufferNav_->get_data_ptr(), sizeof(std::complex<float>)*bufferNav_->get_number_of_elements());

		// kSpace
		GadgetContainerMessage<hoNDArray<std::complex<float>>>* tmp_m3 = new GadgetContainerMessage<hoNDArray<std::complex<float>>>();

		// concatenate data
		tmp_m2->cont(tmp_m3);
		
		// create output
		try{ tmp_m3->getObjectPtr()->create(bufferkSpace_->get_dimensions());}
		catch (std::runtime_error &err){
			GADGET_DEBUG_EXCEPTION(err, "CS_Retro: Unable to allocate new image array\m");
			tmp_m1->release();
			return -1;
		}

		// copy data
		memcpy(tmp_m3->getObjectPtr()->get_data_ptr(), bufferkSpace_->get_data_ptr(), sizeof(std::complex<float>)*bufferkSpace_->get_number_of_elements());

		// put on stream
		if (this->next()->putq(tmp_m1) < 0) {
    		return GADGET_FAIL;
		}
		GADGET_DEBUG2("global PE: %i, PA: %i\n", CS_GlobalVar::instance()->vPE_.size(), CS_GlobalVar::instance()->vPA_.size());

		// save indices to disk
		hoNDArray<float> tmp(CS_GlobalVar::instance()->vNavInd_.size());
		for (int iI = 0; iI < CS_GlobalVar::instance()->vNavInd_.size(); iI++){
			tmp(iI) = CS_GlobalVar::instance()->vNavInd_.at(iI);
		}
		//save_array(tmp, "D:/Martin/aInd");

		return GADGET_OK;
	}
}

// copy all AcquisitionHeader values - from "m1" to "tmp_m1"
int CS_Retro_AccumulatorGadget::copy_header(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *tmp_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1){
	tmp_m1->getObjectPtr()->acquisition_time_stamp		= m1->getObjectPtr()->acquisition_time_stamp;
	tmp_m1->getObjectPtr()->active_channels				= m1->getObjectPtr()->active_channels;
	tmp_m1->getObjectPtr()->available_channels			= m1->getObjectPtr()->available_channels;
	tmp_m1->getObjectPtr()->center_sample				= m1->getObjectPtr()->center_sample;
	tmp_m1->getObjectPtr()->channel_mask[0]				= m1->getObjectPtr()->channel_mask[0];
	tmp_m1->getObjectPtr()->channel_mask[1]				= m1->getObjectPtr()->channel_mask[1];
	tmp_m1->getObjectPtr()->channel_mask[2]				= m1->getObjectPtr()->channel_mask[2];
	tmp_m1->getObjectPtr()->channel_mask[3]				= m1->getObjectPtr()->channel_mask[3];
	tmp_m1->getObjectPtr()->channel_mask[4]				= m1->getObjectPtr()->channel_mask[4];
	tmp_m1->getObjectPtr()->channel_mask[5]				= m1->getObjectPtr()->channel_mask[5];
	tmp_m1->getObjectPtr()->channel_mask[6]				= m1->getObjectPtr()->channel_mask[6];
	tmp_m1->getObjectPtr()->channel_mask[7]				= m1->getObjectPtr()->channel_mask[7];
	tmp_m1->getObjectPtr()->channel_mask[8]				= m1->getObjectPtr()->channel_mask[8];
	tmp_m1->getObjectPtr()->channel_mask[9]				= m1->getObjectPtr()->channel_mask[9];
	tmp_m1->getObjectPtr()->channel_mask[10]			= m1->getObjectPtr()->channel_mask[10];
	tmp_m1->getObjectPtr()->channel_mask[11]			= m1->getObjectPtr()->channel_mask[11];
	tmp_m1->getObjectPtr()->channel_mask[12]			= m1->getObjectPtr()->channel_mask[12];
	tmp_m1->getObjectPtr()->channel_mask[13]			= m1->getObjectPtr()->channel_mask[13];
	tmp_m1->getObjectPtr()->channel_mask[14]			= m1->getObjectPtr()->channel_mask[14];
	tmp_m1->getObjectPtr()->channel_mask[15]			= m1->getObjectPtr()->channel_mask[15];
	tmp_m1->getObjectPtr()->discard_post				= m1->getObjectPtr()->discard_post;
	tmp_m1->getObjectPtr()->discard_pre					= m1->getObjectPtr()->discard_pre;
	tmp_m1->getObjectPtr()->encoding_space_ref			= m1->getObjectPtr()->encoding_space_ref;
	tmp_m1->getObjectPtr()->flags						= m1->getObjectPtr()->flags;
	tmp_m1->getObjectPtr()->idx.average					= m1->getObjectPtr()->idx.average;
	tmp_m1->getObjectPtr()->idx.contrast				= m1->getObjectPtr()->idx.contrast;
	tmp_m1->getObjectPtr()->idx.kspace_encode_step_1	= m1->getObjectPtr()->idx.kspace_encode_step_1;
	tmp_m1->getObjectPtr()->idx.kspace_encode_step_2	= m1->getObjectPtr()->idx.kspace_encode_step_2;
	tmp_m1->getObjectPtr()->idx.phase					= m1->getObjectPtr()->idx.phase;
	tmp_m1->getObjectPtr()->idx.repetition				= m1->getObjectPtr()->idx.repetition;
	tmp_m1->getObjectPtr()->idx.segment					= m1->getObjectPtr()->idx.segment;
	tmp_m1->getObjectPtr()->idx.set						= m1->getObjectPtr()->idx.set;
	tmp_m1->getObjectPtr()->idx.slice					= m1->getObjectPtr()->idx.slice;
	tmp_m1->getObjectPtr()->idx.user[0]					= m1->getObjectPtr()->idx.user[0];
	tmp_m1->getObjectPtr()->idx.user[1]					= m1->getObjectPtr()->idx.user[1];
	tmp_m1->getObjectPtr()->idx.user[2]					= m1->getObjectPtr()->idx.user[2];
	tmp_m1->getObjectPtr()->idx.user[3]					= m1->getObjectPtr()->idx.user[3];
	tmp_m1->getObjectPtr()->idx.user[4]					= m1->getObjectPtr()->idx.user[4];
	tmp_m1->getObjectPtr()->idx.user[5]					= m1->getObjectPtr()->idx.user[5];
	tmp_m1->getObjectPtr()->idx.user[6]					= m1->getObjectPtr()->idx.user[6];
	tmp_m1->getObjectPtr()->idx.user[7]					= m1->getObjectPtr()->idx.user[7];
	tmp_m1->getObjectPtr()->measurement_uid				= m1->getObjectPtr()->measurement_uid;
	tmp_m1->getObjectPtr()->number_of_samples			= m1->getObjectPtr()->number_of_samples;
	tmp_m1->getObjectPtr()->patient_table_position[0]	= m1->getObjectPtr()->patient_table_position[0];
	tmp_m1->getObjectPtr()->patient_table_position[1]	= m1->getObjectPtr()->patient_table_position[1];
	tmp_m1->getObjectPtr()->patient_table_position[2]	= m1->getObjectPtr()->patient_table_position[2];
	tmp_m1->getObjectPtr()->phase_dir[0]				= m1->getObjectPtr()->phase_dir[0];
	tmp_m1->getObjectPtr()->phase_dir[1]				= m1->getObjectPtr()->phase_dir[1];
	tmp_m1->getObjectPtr()->phase_dir[2]				= m1->getObjectPtr()->phase_dir[2];
	tmp_m1->getObjectPtr()->physiology_time_stamp[0]	= m1->getObjectPtr()->physiology_time_stamp[0];
	tmp_m1->getObjectPtr()->physiology_time_stamp[1]	= m1->getObjectPtr()->physiology_time_stamp[1];
	tmp_m1->getObjectPtr()->physiology_time_stamp[2]	= m1->getObjectPtr()->physiology_time_stamp[2];
	tmp_m1->getObjectPtr()->physiology_time_stamp[3]	= m1->getObjectPtr()->physiology_time_stamp[3];
	tmp_m1->getObjectPtr()->physiology_time_stamp[4]	= m1->getObjectPtr()->physiology_time_stamp[4];
	tmp_m1->getObjectPtr()->physiology_time_stamp[5]	= m1->getObjectPtr()->physiology_time_stamp[5];
	tmp_m1->getObjectPtr()->physiology_time_stamp[6]	= m1->getObjectPtr()->physiology_time_stamp[6];
	tmp_m1->getObjectPtr()->physiology_time_stamp[7]	= m1->getObjectPtr()->physiology_time_stamp[7];
	tmp_m1->getObjectPtr()->position[0]					= m1->getObjectPtr()->position[0];
	tmp_m1->getObjectPtr()->position[1]					= m1->getObjectPtr()->position[1];
	tmp_m1->getObjectPtr()->position[2]					= m1->getObjectPtr()->position[2];
	tmp_m1->getObjectPtr()->read_dir[0]					= m1->getObjectPtr()->read_dir[0];
	tmp_m1->getObjectPtr()->read_dir[1]					= m1->getObjectPtr()->read_dir[1];
	tmp_m1->getObjectPtr()->read_dir[2]					= m1->getObjectPtr()->read_dir[2];
	tmp_m1->getObjectPtr()->sample_time_us				= m1->getObjectPtr()->sample_time_us;
	tmp_m1->getObjectPtr()->scan_counter				= m1->getObjectPtr()->scan_counter;
	tmp_m1->getObjectPtr()->slice_dir[0]				= m1->getObjectPtr()->slice_dir[0];
	tmp_m1->getObjectPtr()->slice_dir[1]				= m1->getObjectPtr()->slice_dir[1];
	tmp_m1->getObjectPtr()->slice_dir[2]				= m1->getObjectPtr()->slice_dir[2];
	tmp_m1->getObjectPtr()->trajectory_dimensions		= m1->getObjectPtr()->trajectory_dimensions;
	tmp_m1->getObjectPtr()->user_float[0]				= m1->getObjectPtr()->user_float[0];
	tmp_m1->getObjectPtr()->user_float[1]				= m1->getObjectPtr()->user_float[1];
	tmp_m1->getObjectPtr()->user_float[2]				= m1->getObjectPtr()->user_float[2];
	tmp_m1->getObjectPtr()->user_float[3]				= m1->getObjectPtr()->user_float[3];
	tmp_m1->getObjectPtr()->user_float[4]				= m1->getObjectPtr()->user_float[4];
	tmp_m1->getObjectPtr()->user_float[5]				= m1->getObjectPtr()->user_float[5];
	tmp_m1->getObjectPtr()->user_float[6]				= m1->getObjectPtr()->user_float[6];
	tmp_m1->getObjectPtr()->user_float[7]				= m1->getObjectPtr()->user_float[7];
	tmp_m1->getObjectPtr()->user_int[0]					= m1->getObjectPtr()->user_int[0];
	tmp_m1->getObjectPtr()->user_int[1]					= m1->getObjectPtr()->user_int[1];
	tmp_m1->getObjectPtr()->user_int[2]					= m1->getObjectPtr()->user_int[2];
	tmp_m1->getObjectPtr()->user_int[3]					= m1->getObjectPtr()->user_int[3];
	tmp_m1->getObjectPtr()->user_int[4]					= m1->getObjectPtr()->user_int[4];
	tmp_m1->getObjectPtr()->user_int[5]					= m1->getObjectPtr()->user_int[5];
	tmp_m1->getObjectPtr()->user_int[6]					= m1->getObjectPtr()->user_int[6];
	tmp_m1->getObjectPtr()->user_int[7]					= m1->getObjectPtr()->user_int[7];
	tmp_m1->getObjectPtr()->version						= m1->getObjectPtr()->version;
	
	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_AccumulatorGadget)

}

/*
		//------------------------------------------------------------------------------------
		//------------------ get navigator signal from navigator array -----------------------
		//------------------------------------------------------------------------------------
		if (getNav2D(*bufferNav_)){
			GADGET_DEBUG1("CS_Retro: Error occurred in function getNav2D(...)\n");
			return GADGET_FAIL;
		}

		//------------------------------------------------------------------------------------
		//----- reorder k-space depending on the number of phases and navigator signal -------
		//------------------------------------------------------------------------------------
		hoNDArray<std::complex<float>> kSpaceOrdered;
		if (!reorder_kSpace(*bufferkSpace_, kSpaceOrdered)){
			GADGET_DEBUG1("CS_Retro: Error occurred in function reorder_kSpace(...)\n");
			return GADGET_FAIL;
		}

		//------------------------------------------------------------------------------------
		// copy ordered kSpace to a new array, create new header information and put on stream
		//------------------------------------------------------------------------------------
		// create new image header
		GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();

		// initialize the image header
		memset(cm1->getObjectPtr(),0,sizeof(ISMRMRD::ImageHeader));
    
		// initialize flags
		cm1->getObjectPtr()->flags = 0;
	
		// create new GadgetContainer
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 = new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
    
		// concatenate data with header
		cm1->cont(cm2);
	
		// create output
		try{ cm2->getObjectPtr()->create(kSpaceOrdered.get_dimensions());}
		catch (std::runtime_error &err){
			GADGET_DEBUG_EXCEPTION(err, "CS_Retro: Unable to allocate new image array\m");
			cm1->release();
			return -1;
		}

		// set 4D flag
		cm1->getObjectPtr()->user_int[0] = 7;

		// copy data
		memcpy(cm2->getObjectPtr()->get_data_ptr(), kSpaceOrdered.get_data_ptr(), sizeof(std::complex<float>)*kSpaceOrdered.get_number_of_elements());

		// set header information
		cm1->getObjectPtr()->user_int[0] = 7;
		cm1->getObjectPtr()->matrix_size[0]     = (uint16_t)kSpaceOrdered.get_size(0);
		cm1->getObjectPtr()->matrix_size[1]     = (uint16_t)kSpaceOrdered.get_size(1);
		cm1->getObjectPtr()->matrix_size[2]     = (uint16_t)kSpaceOrdered.get_size(2);
		cm1->getObjectPtr()->field_of_view[0]   = field_of_view_[0];
		cm1->getObjectPtr()->field_of_view[1]   = field_of_view_[1];
		cm1->getObjectPtr()->field_of_view[2]   = field_of_view_[2];
		cm1->getObjectPtr()->channels           = (uint16_t)m1->getObjectPtr()->active_channels;
		cm1->getObjectPtr()->slice				= m1->getObjectPtr()->idx.slice;
		memcpy(cm1->getObjectPtr()->position,m1->getObjectPtr()->position,sizeof(float)*3);
		memcpy(cm1->getObjectPtr()->read_dir,m1->getObjectPtr()->read_dir,sizeof(float)*3);
		memcpy(cm1->getObjectPtr()->phase_dir,m1->getObjectPtr()->phase_dir,sizeof(float)*3);
		memcpy(cm1->getObjectPtr()->slice_dir,m1->getObjectPtr()->slice_dir, sizeof(float)*3);
		memcpy(cm1->getObjectPtr()->patient_table_position,m1->getObjectPtr()->patient_table_position, sizeof(float)*3);
		cm1->getObjectPtr()->image_data_type = ISMRMRD::DATA_COMPLEX_FLOAT;
		cm1->getObjectPtr()->image_index = (uint16_t)(++image_counter_);
		cm1->getObjectPtr()->image_series_index = (uint16_t)image_series_;

		// set user values, if Compressed Sensing is active
		if(this->get_bool_value("CS_on") == true){
			cm1->getObjectPtr()->user_float[0] = fCSAcc_;
			cm1->getObjectPtr()->user_float[1] = fFullySa_/100;
			cm1->getObjectPtr()->user_float[2] = fPartialFourierVal_;
			cm1->getObjectPtr()->user_float[3] = fLQ_;
			cm1->getObjectPtr()->user_float[4] = fLESPReSSo_;

			cm1->getObjectPtr()->user_int[1] = iBodyRegion_;
			cm1->getObjectPtr()->user_int[2] = iSamplingType_;
			cm1->getObjectPtr()->user_int[3] = iVDMap_;
			cm1->getObjectPtr()->user_int[4] = iESPReSSoDirection_;
		}

		// correct readout, phase encoding, slice encoding direction
		if (m1->getObjectPtr()->read_dir[0] < 0 || m1->getObjectPtr()->read_dir[1] < 0 || m1->getObjectPtr()->read_dir[2] < 0){
			GADGET_DEBUG1("readout direction is inverse - flip data before start streaming\n");
			flip_array(*cm2->getObjectPtr(), 0);
		}
		if (m1->getObjectPtr()->phase_dir[0] < 0 || m1->getObjectPtr()->phase_dir[1] < 0 || m1->getObjectPtr()->phase_dir[2] < 0){
			GADGET_DEBUG1("phase encoding direction is inverse - flip data before start streaming\n");
			flip_array(*cm2->getObjectPtr(), 1);
		}
		if (m1->getObjectPtr()->slice_dir[0] < 0 || m1->getObjectPtr()->slice_dir[1] < 0 || m1->getObjectPtr()->slice_dir[2] < 0){
			GADGET_DEBUG1("partition encoding direction is inverse - flip data before start streaming\n");
			flip_array(*cm2->getObjectPtr(), 2);
		}

		// put on stream
		if (this->next()->putq(cm1) < 0) {
    		return GADGET_FAIL;
		}
		*/
