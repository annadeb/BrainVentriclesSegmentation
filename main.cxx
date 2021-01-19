#include <iostream> // std::cout, std::cin
#include <string> // std::string
#include <itkImage.h> // itk::Image
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageSeriesReader.h>
#include <itkImageSeriesWriter.h>
#include <itkNumericSeriesFileNames.h>
#include <itkGDCMSeriesFileNames.h> 
#include <itkGDCMImageIO.h>
#include<itkGDCMImageIO.h>
#include<itkVTKImageIO.h>
#include<itkIntensityWindowingImageFilter.h>
#include<itkThresholdImageFilter.h>
#include<itkBinaryThresholdImageFilter.h>
#include<itkOtsuMultipleThresholdsImageFilter.h>
#include<itkBinaryBallStructuringElement.h>
#include<itkErodeObjectMorphologyImageFilter.h>
#include<itkBinaryErodeImageFilter.h>
#include<itkBinaryDilateImageFilter.h>
#include<itkConnectedComponentImageFilter.h>
#include<itkGradientMagnitudeImageFilter.h>
#include <itkCropImageFilter.h>
#include<itkConnectedThresholdImageFilter.h>
#include<itkBinaryImageToLabelMapFilter.h>
#include<itkLabelMapToLabelImageFilter.h>
#include<itkLabelStatisticsImageFilter.h>
#include<itkLabelShapeKeepNObjectsImageFilter.h>
#include<itkRescaleIntensityImageFilter.h>
#include<itkBinaryFillholeImageFilter.h>
#include<itkVotingBinaryIterativeHoleFillingImageFilter.h>
#include<itkMultiplyImageFilter.h>
#include<itkMedianImageFilter.h>
#include<itkMinimumMaximumImageCalculator.h>
#include<itkBinaryMorphologicalOpeningImageFilter.h>
#include<itkConfidenceConnectedImageFilter.h>
//#include<itkCurvatureFlowImageFilter.h>
#include<itkNeighborhoodConnectedImageFilter.h>
//#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkInvertIntensityImageFilter.h>
#include<itkRelabelComponentImageFilter.h>
#include<itkFlipImageFilter.h>
#include<itkAddImageFilter.h>
#include<itkJoinImageFilter.h>
#include "itkLabelGeometryImageFilter.h"





using PixelType = signed short;
using ImageType = itk::Image<PixelType, 2>;
using BallType = itk::BinaryBallStructuringElement<short, 3>;
//using InputPixelType = float;
//using OutputPixelType = int;

using Image3DType = itk::Image<PixelType, 3>;
using ImageIOType = itk::GDCMImageIO;
//using InputImageType = itk::Image< InputPixelType, 3 >;
//using OutputImageType = itk::Image< OutputPixelType, 3 >;
using SeriesReaderType = itk::ImageSeriesReader<Image3DType>;
using Series2DWriterType = itk::ImageSeriesWriter<Image3DType, ImageType>;
using Series3DWriterType = itk::ImageSeriesWriter<Image3DType, Image3DType>;

using Writer3Dtype = itk::ImageFileWriter<Image3DType>;
using ReaderType = itk::ImageFileReader<ImageType>; 
using Reader3DType = itk::ImageFileReader<Image3DType>;

using WriterType = itk::ImageFileWriter<ImageType>;

using ConnectedComponent = itk::ConnectedComponentImageFilter<Image3DType, Image3DType, Image3DType>;
using NumSeriesFileNames = itk::NumericSeriesFileNames;

using ConnectedThreshold = itk::ConnectedThresholdImageFilter<Image3DType, Image3DType>;

int main() // glowna funkcja programu
{
	try {
		ReaderType::Pointer reader = ReaderType::New();
		Reader3DType::Pointer reader3D = Reader3DType::New();
		Reader3DType::Pointer reader3Db = Reader3DType::New();
		Reader3DType::Pointer reader3Dc = Reader3DType::New();
		Reader3DType::Pointer reader3Dd = Reader3DType::New();
		WriterType::Pointer writer = WriterType::New();
		ImageType::Pointer image = ImageType::New();

		SeriesReaderType::Pointer seriesReader = SeriesReaderType::New();
		Series2DWriterType::Pointer series2DWriter = Series2DWriterType::New();
		itk::GDCMSeriesFileNames::Pointer namesGen = itk::GDCMSeriesFileNames::New();
		NumSeriesFileNames::Pointer numSeriesFileNames = NumSeriesFileNames::New();
		Image3DType::Pointer image3D = Image3DType::New();
		Image3DType::Pointer image3D_bin = Image3DType::New();
		Image3DType::Pointer image3D_mnozenie = Image3DType::New();

		ImageIOType::Pointer dicomIO = ImageIOType::New();
		Writer3Dtype::Pointer writer3D = Writer3Dtype::New();
		Series3DWriterType::Pointer series3DWriter = Series3DWriterType::New();
		ConnectedComponent::Pointer CCImageFilter = ConnectedComponent::New();
		itk::GDCMImageIO::Pointer gdcmImageIO = itk::GDCMImageIO::New();


		//namesGen->SetDirectory("../dane_glowa/glowa/mozg md/Head_Neck_Standard - 1/t2_tse_tra_6");
		//namesGen->SetDirectory("../dane_glowa/glowa/mozg sd/Head_Neck_Standard - 1/t2_tse_tra_6");
		
		namesGen->SetDirectory("../dane_glowa/glowa/mozg zb/Head_Neck_Standard - 153275/t2_tse_tra_5");
		itk::GDCMSeriesFileNames::SeriesUIDContainerType seriesUIDs = namesGen->GetSeriesUIDs();
		itk::GDCMSeriesFileNames::FileNamesContainerType fileNames = namesGen->GetFileNames(seriesUIDs[0]);
		for (size_t i = 0; i < seriesUIDs.size(); i++) {
			std::cout << seriesUIDs[i] << std::endl;
		}
		//Czytanie serii obrazowej 
		seriesReader->SetFileNames(namesGen->GetFileNames(seriesUIDs[0]));
		seriesReader->SetImageIO(gdcmImageIO);
		image3D = seriesReader->GetOutput();
		seriesReader->Update();


		//std::cout << image3D;
		//std::cout << seriesReader->GetMetaDataDictionaryArray();

		//Zapisywanie serii obrazowej 

		//numSeriesFileNames->SetSeriesFormat("..\\wyniki\\IMG\%05d.dcm");
		//numSeriesFileNames->SetStartIndex(1);
		//numSeriesFileNames->SetEndIndex(image3D->GetLargestPossibleRegion().GetSize()[2]);

		//series2DWriter->SetFileNames(numSeriesFileNames->GetFileNames());
	 //   series2DWriter->SetImageIO(gdcmImageIO);
		////series2DWriter->SetMetaDataDictionaryArray(seriesReader->GetMetaDataDictionaryArray());
		//series2DWriter->SetInput(image3D);
		//series2DWriter->Update();
		//
		series3DWriter->SetInput(image3D);
		series3DWriter->SetFileName("..\\wyniki\\img3D.vtk");
		series3DWriter->Update();


		//binaryzacja
		using FilterType = itk::BinaryThresholdImageFilter<Image3DType, Image3DType>;
		FilterType::Pointer thresholder = FilterType::New();
		thresholder->SetInput(image3D);
		thresholder->SetInsideValue(1);
		thresholder->SetOutsideValue(0);
		thresholder->SetLowerThreshold(300);

		series3DWriter->SetInput(thresholder->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\img3D_bin_calosc.vtk");
		series3DWriter->Update();
		image3D_bin = thresholder->GetOutput();
		//zalewanie obszaru
		/*using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<Image3DType, Image3DType>;

		ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
		connected->SetInput(image3D_bin);
		connected->SetBackgroundValue(1);
			connected->Update();

		std::cout << "Number of objects: " << connected->GetObjectCount() << std::endl;
		series3DWriter->SetInput(connected->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\img3D_zalane.vtk");
		series3DWriter->Update();*/

		//dylatacja
		BallType::SizeType rad;
		rad[0] = 8;
		rad[1] = 8;
		rad[2] = 8;
		using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, Image3DType::ImageDimension>;
		StructuringElementType structuringElement;
		structuringElement.SetRadius(rad);
		structuringElement.CreateStructuringElement();

		using FilterDilateType = itk::BinaryDilateImageFilter<Image3DType, Image3DType, StructuringElementType>;
		FilterDilateType::Pointer dilateFilter = FilterDilateType::New();
		dilateFilter->SetInput(image3D_bin);
		dilateFilter->SetKernel(structuringElement);
		dilateFilter->SetBackgroundValue(0);
		dilateFilter->SetForegroundValue(1);
		dilateFilter->Update();

		series3DWriter->SetInput(dilateFilter->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\img3D_dilate.vtk");
		series3DWriter->Update();

		//erozja 
		rad[0] = 38;
		rad[1] = 38;
		rad[2] = 38;
		structuringElement.SetRadius(rad);
		structuringElement.CreateStructuringElement();

		using FilterErodeType = itk::BinaryErodeImageFilter<Image3DType, Image3DType, StructuringElementType>;
		FilterErodeType::Pointer erodeFilter = FilterErodeType::New();
		erodeFilter->SetInput(dilateFilter->GetOutput());
		erodeFilter->SetKernel(structuringElement);
		erodeFilter->SetBackgroundValue(0);
		erodeFilter->SetForegroundValue(1);
		erodeFilter->Update();

		series3DWriter->SetInput(erodeFilter->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\img3D_erozja.vtk");
		series3DWriter->Update();

		//mno¿enie
		//image3D* erodeFilter->GetOutput();
		using MultiplyType = itk::MultiplyImageFilter< Image3DType, Image3DType, Image3DType>;
		MultiplyType::Pointer multiply = MultiplyType::New();
		multiply->SetInput1(image3D);
		multiply->SetInput2(erodeFilter->GetOutput());
		image3D_mnozenie = multiply->GetOutput();;
		multiply->Update();
		series3DWriter->SetInput(image3D_mnozenie);
		series3DWriter->SetFileName("..\\wyniki\\img3D_mnozenie.vtk");
		series3DWriter->Update();

		//filter medianowy 
		//using MedianFilterType = itk::MedianImageFilter<Image3DType, Image3DType>;
		//MedianFilterType::Pointer medianFilter = MedianFilterType::New();
		//MedianFilterType::InputSizeType radius;
		//radius.Fill(1);
		//medianFilter->SetRadius(radius);
		//medianFilter->SetInput(multiply->GetOutput());	//using FillholeFilterType = itk::BinaryFillholeImageFilter<Image3DType>;
		//series3DWriter->SetInput(medianFilter->GetOutput());
		//series3DWriter->SetFileName("..\\wyniki\\img3D_filtr_medianowy.vtk");
		//series3DWriter->Update();

		//binaryzacja -> poszukiwanie max
		thresholder->SetInput(multiply->GetOutput());
		thresholder->SetInsideValue(1);
		thresholder->SetOutsideValue(0);
		thresholder->SetLowerThreshold(600);
		thresholder->SetUpperThreshold(800);
		series3DWriter->SetInput(thresholder->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\img3D_bin_komory.vtk");
		series3DWriter->Update();

		//otwarcie
		using OpeningFilterType = itk::BinaryMorphologicalOpeningImageFilter<Image3DType, Image3DType, StructuringElementType>;
		int radiusOpen = 1;
		structuringElement.SetRadius(radiusOpen);
		structuringElement.CreateStructuringElement();

		OpeningFilterType::Pointer openFilter = OpeningFilterType::New();
		openFilter->SetInput(thresholder->GetOutput());
		openFilter->SetKernel(structuringElement);
		openFilter->SetBackgroundValue(0);
		openFilter->SetForegroundValue(1);
		openFilter->Update();

		series3DWriter->SetInput(openFilter->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\img3D_otwarcie.vtk");
		series3DWriter->Update();

		//pierwszy niezerowy punkt obrazu od czubka g³owy w dó³
		//TODO

		//rozrost
	/*	using CurvatureFlowImageFilterType =itk::CurvatureFlowImageFilter< Image3DType, Image3DType >;
		CurvatureFlowImageFilterType::Pointer smoothing =CurvatureFlowImageFilterType::New();
		smoothing->SetInput(image3D);
		smoothing->SetNumberOfIterations(5);
		smoothing->SetTimeStep(0.125);*/

		using ConnCompFilterType = itk::ConnectedComponentImageFilter<Image3DType, Image3DType>;

		ConnCompFilterType::Pointer connComp1 = ConnCompFilterType::New();
		connComp1->SetInput(openFilter->GetOutput());
		connComp1->SetBackgroundValue(0);
		connComp1->Update();
		series3DWriter->SetInput(connComp1->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\img3D_etykiety.vtk");
		series3DWriter->Update();


		using LabelGeometryImageFilterType = itk::LabelGeometryImageFilter<Image3DType>;
		LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
		labelGeometryImageFilter->SetInput(connComp1->GetOutput());
		//labelGeometryImageFilter->SetIntensityInput(connComp->GetOutput());

		// These generate optional outputs.
		labelGeometryImageFilter->CalculatePixelIndicesOn();
		labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();

		labelGeometryImageFilter->Update();
		LabelGeometryImageFilterType::LabelsType labelsAll= labelGeometryImageFilter->GetLabels();
		
		std::cout << labelsAll.size() << std::endl;
		std::cout << labelGeometryImageFilter->GetCentroid(labelsAll.size()-1) << std::endl;
		
		LabelGeometryImageFilterType::LabelPointType startPoint = labelGeometryImageFilter->GetCentroid(labelsAll.size() - 1);
		//LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
		
		std::cout << startPoint[0] << std::endl;
		std::cout << startPoint[1] << std::endl;
		std::cout << startPoint[2] << std::endl;
	//===============================================================================================================

		using ConnectedFilterType =itk::ConfidenceConnectedImageFilter<Image3DType, Image3DType>;
		ConnectedFilterType::Pointer confidenceConnected = ConnectedFilterType::New();
		confidenceConnected->SetInput(image3D); // daæ tu obraz po mno¿eniu zeby bylo przyciête
		//confidenceConnected->SetMultiplier(0.57);  //Multiplier = 0.57, iteracje 1
		//confidenceConnected->SetNumberOfIterations(1);
		confidenceConnected->SetMultiplier(1.3);  //Multiplier = 0.57, iteracje 1
		confidenceConnected->SetNumberOfIterations(1);
		confidenceConnected->SetReplaceValue(99);
		confidenceConnected->SetInitialNeighborhoodRadius(2);
		ConnectedFilterType::IndexType index;
		index[0] = startPoint[0]+10; // 152;//129;//152;
		index[1] = startPoint[1]; //168;//129;// 168;
		index[2] = startPoint[2];//48;// 48;

		//index[0] = 111;//117;//129;//152;
		//index[1] = 130;// 158;//129;// 168;
		//index[2] = 47;// 48;
		//ConnectedFilterType::IndexType index2;
		//index2[0] = 111;//129;//152;
		//index2[1] = 168;//129;// 168;
		//index2[2] = 48;// 48;
	confidenceConnected->AddSeed(index);
	//confidenceConnected->AddSeed(index2);
	confidenceConnected->Update();
	series3DWriter->SetInput(confidenceConnected->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_rozrost-thresh.vtk");
	series3DWriter->Update();
	auto pixelVal=image3D->GetPixel(index);
	
	std::cout << image3D_mnozenie->GetPixel(index) << std::endl;

	//multiply->SetInput1(confidenceConnected->GetOutput());
	//multiply->SetInput2(image3D_mnozenie);
	////image3D_mnozenie = multiply->GetOutput();;
	//multiply->Update();
	//series3DWriter->SetInput(multiply->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_mnozenie_v2.vtk");
	//series3DWriter->Update();


	using InvertIntensityImageFilterType = itk::InvertIntensityImageFilter<Image3DType>;

	InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
	invertIntensityFilter->SetInput(confidenceConnected->GetOutput());
	invertIntensityFilter->SetMaximum(1);
	invertIntensityFilter->Update();
	series3DWriter->SetInput(invertIntensityFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_invert.vtk");
	series3DWriter->Update();

	using WindowingImageFilter = itk::IntensityWindowingImageFilter<Image3DType>;
	WindowingImageFilter::Pointer windowingImageFilter = WindowingImageFilter::New();
	windowingImageFilter->SetInput(invertIntensityFilter->GetOutput());
	windowingImageFilter->SetWindowMaximum(1);
	windowingImageFilter->SetWindowMinimum(-1);
	windowingImageFilter->SetOutputMinimum(0);
	windowingImageFilter->SetOutputMaximum(1);
	windowingImageFilter->Update();
	series3DWriter->SetInput(windowingImageFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_window.vtk");
	series3DWriter->Update();

	

	Image3DType::Pointer windowImage = Image3DType::New();
	Image3DType::Pointer erode = Image3DType::New();

	

	reader3D->SetFileName("..\\wyniki\\img3D_window.vtk");
	reader3D->Update();
	windowImage = reader3D->GetOutput();

	reader3Db->SetFileName("..\\wyniki\\img3D_erozja.vtk");
	reader3Db->Update();
	erode = reader3Db->GetOutput();

	index[0] = 111;//117;//129;//152;
	index[1] = 130;// 158;//129;// 168;
	index[2] = 47;// 48;

	std::cout << windowingImageFilter->GetOutput()->GetPixel(index) << std::endl;
	std::cout << reader3Db->GetOutput()->GetPixel(index) << std::endl;

	multiply->SetInput1(windowImage);
	multiply->SetInput2(erode);
	//image3D_mnozenie = multiply->GetOutput();;
	multiply->Update();
	series3DWriter->SetInput(multiply->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_mnozenie_komory.vtk");
	series3DWriter->Update();





	using ConnCompFilterType = itk::ConnectedComponentImageFilter<Image3DType, Image3DType>;

	ConnCompFilterType::Pointer connComp = ConnCompFilterType::New();
	connComp->SetInput(multiply->GetOutput());
	connComp->SetBackgroundValue(0);
	connComp->Update();

	series3DWriter->SetInput(connComp->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_label.vtk");
	series3DWriter->Update();


	using RelabelComponentFilterType = itk::RelabelComponentImageFilter<Image3DType, Image3DType>;
	RelabelComponentFilterType::Pointer relabel = RelabelComponentFilterType::New();
	relabel->SetInput(connComp->GetOutput());
	relabel->SetMinimumObjectSize(150);
	
	relabel->Update();

	series3DWriter->SetInput(relabel->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_relabel.vtk");
	series3DWriter->Update();


	using FilterType = itk::BinaryThresholdImageFilter<Image3DType, Image3DType>;
	FilterType::Pointer thresholder2 = FilterType::New();
	thresholder2->SetInput(relabel->GetOutput());
	thresholder2->SetLowerThreshold(2);
	thresholder2->SetUpperThreshold(2);
	thresholder2->Update();

	series3DWriter->SetInput(thresholder2->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_binaryzacja1komora.vtk");
	series3DWriter->Update();
	
	//flip
	using FlipImageFilterType = itk::FlipImageFilter<Image3DType>;

	FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
	flipFilter->SetInput(thresholder2->GetOutput());


	FlipImageFilterType::FlipAxesArrayType flipAxes;

	flipAxes[0] = true;
	flipAxes[1] = false;
	flipAxes[2] = false;
	flipFilter->SetFlipAxes(flipAxes);
	flipFilter->Update();

	series3DWriter->SetInput(flipFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\img3D_flip.vtk");
	series3DWriter->Update();

	//add
	//reader3Dc->SetFileName("..\\wyniki\\img3D_binaryzacja1komora.vtk");
	////reader3Dc->SetFileName("..\\wyniki\\img3D.vtk");
	//reader3Dc->Update();
	//Image3DType::Pointer komora13D = Image3DType::New();
	//komora13D = reader3Dc->GetOutput();

	//reader3Dd->SetFileName("..\\wyniki\\img3D_flip.vtk");
	//reader3Dd->Update();
	//Image3DType::Pointer flip = Image3DType::New();
	//flip = reader3Dd->GetOutput();

	//using AddImageFilterType = itk::AddImageFilter<Image3DType, Image3DType>;

	//AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

	///*using JoinImageFilterType = itk::JoinImageFilter<Image3DType, Image3DType>;

	//JoinImageFilterType::Pointer addFilter = JoinImageFilterType::New();*/
	//addFilter->SetInput1(flip);
	//addFilter->SetInput2(komora13D);
	//

	//series3DWriter->SetInput(addFilter->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_dodawanie.vtk");
	//series3DWriter->Update();

	
	//using AbsoluteValueDifferenceImageFilterType =	itk::AbsoluteValueDifferenceImageFilter<Image3DType, Image3DType, Image3DType>;

	//AbsoluteValueDifferenceImageFilterType::Pointer absoluteValueDifferenceFilter =
	//	AbsoluteValueDifferenceImageFilterType::New();
	//absoluteValueDifferenceFilter->SetInput1(image3D_mnozenie);
	//absoluteValueDifferenceFilter->SetInput2(confidenceConnected->GetOutput());
	//absoluteValueDifferenceFilter->Update();

	//series3DWriter->SetInput(absoluteValueDifferenceFilter->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_difference.vtk");
	//series3DWriter->Update();


	//	using ConnectedFilterType =	itk::ConnectedThresholdImageFilter< Image3DType,Image3DType >;
	//	ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
	//	connectedThreshold->SetInput(image3D);
	//	connectedThreshold->SetLower(500);
	//	connectedThreshold->SetConnectivity(ConnectedFilterType::FullConnectivity);
	//	connectedThreshold->SetUpper(820);
	//	connectedThreshold->SetReplaceValue(255);
	//			ConnectedFilterType::IndexType index;
	//	index[0] = 129;//129;//152;
	//	index[1] = 129;//129;// 168;
	//	index[2] = 48;// 48;
	//	connectedThreshold->SetSeed(index);
	//	connectedThreshold->SetSeed(index);
	//series3DWriter->SetInput(connectedThreshold->GetOutput());
	//series3DWriter->SetFileName("..\\wyniki\\img3D_rozrost-thresh.vtk");
	//series3DWriter->Update();


		//using ConnectedFilterType =	itk::NeighborhoodConnectedImageFilter<Image3DType, Image3DType >;
		//ConnectedFilterType::Pointer neighborhoodConnected	= ConnectedFilterType::New();
		//neighborhoodConnected->SetInput(image3D);
		//neighborhoodConnected->SetLower(0);
		//neighborhoodConnected->SetUpper(200);

		//Image3DType::SizeType radius3D;

		//radius3D[0] = 5;   // two pixels along X
		//radius3D[1] = 5;   // two pixels along Y
		//radius3D[2] =5;

		//neighborhoodConnected->SetRadius(radius3D);
		//
		//neighborhoodConnected->SetReplaceValue(255);

		//ConnectedFilterType::IndexType index;

		/*Image3DType::RegionType region = image3D->GetLargestPossibleRegion();
		Image3DType::SizeType size = region.GetSize();
		std::cout << size << std::endl;*/

		//index[0] = 152;//129;//152;
		//index[1] = 168;//129;// 168;
		//index[2] = 48;// 48;
		//neighborhoodConnected->SetSeed(index);

		/*confidenceConnected->SetSeed(index);
		confidenceConnected->SetInitialNeighborhoodRadius(2);
		confidenceConnected->Update();*/

		

		//series3DWriter->SetInput(neighborhoodConnected->GetOutput());
		//series3DWriter->SetFileName("..\\wyniki\\img3D_rozrost-thresh.vtk");
		//series3DWriter->Update();



		//	using LabelShapeKeepNObjectsImageFilterType = itk::LabelShapeKeepNObjectsImageFilter<OutputImageType>;
		//	LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter =
		//		LabelShapeKeepNObjectsImageFilterType::New();
		//	labelShapeKeepNObjectsImageFilter->SetInput(connected->GetOutput());
		//	labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0);
		//	labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(1);
		//	labelShapeKeepNObjectsImageFilter->SetAttribute(
		//		LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
		//	
		//	using RescaleFilterType = itk::RescaleIntensityImageFilter<OutputImageType, Image3DType>;
		//	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		//	rescaleFilter->SetOutputMinimum(0);
		//	rescaleFilter->SetOutputMaximum(itk::NumericTraits<PixelType>::max());
		//	rescaleFilter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
		//
		//#ifdef ENABLE_QUICKVIEW
		//	QuickView viewer;
		//	viewer.AddImage(image3D, true, itksys::SystemTools::GetFilenameName(argv[1]));
		//
		//	std::stringstream desc;
		//	desc << "Largest object of " << connected->GetObjectCount() << " objects";
		//	viewer.AddImage(rescaleFilter->GetOutput(), true, desc.str());
		//
		//	viewer.Visualize();
		//#endif
		//
		//	//
		//	using BinaryImageToLabelMapFilterType = itk::BinaryImageToLabelMapFilter<Image3DType>;
		//	BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
		//	binaryImageToLabelMapFilter->SetInput(openFilter->GetOutput());
		//	binaryImageToLabelMapFilter->Update();
		////
		//	using LabelMapToLabelImageFilterType =
		//		itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, Image3DType>;
		//	LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
		//	labelMapToLabelImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
		//	labelMapToLabelImageFilter->Update();
		////
		//	using LabelStatisticsImageFilterType = itk::LabelStatisticsImageFilter<Image3DType, Image3DType>;
		//	LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
		//	labelStatisticsImageFilter->SetLabelInput(labelMapToLabelImageFilter->GetOutput());
		//	labelStatisticsImageFilter->SetInput(openFilter->GetOutput());
		//	labelStatisticsImageFilter->Update();
		//	labelStatisticsImageFilter->GetRegion(labelStatisticsImageFilter->GetValidLabelValues()[0]);

		//
		//	std::cout << "Number of labels: " << labelStatisticsImageFilter->GetNumberOfLabels() << std::endl;
		//	std::cout << std::endl;
		//
		//	using LabelPixelType = LabelStatisticsImageFilterType::LabelPixelType;
		//
		//	for (auto vIt = labelStatisticsImageFilter->GetValidLabelValues().begin();
		//		vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
		//		++vIt)
		//	{
		//		if (labelStatisticsImageFilter->HasLabel(*vIt))
		//		{
		//			LabelPixelType labelValue = *vIt;
		//			std::cout << "min: " << labelStatisticsImageFilter->GetMinimum(labelValue) << std::endl;
		//			std::cout << "max: " << labelStatisticsImageFilter->GetMaximum(labelValue) << std::endl;
		//			std::cout << "median: " << labelStatisticsImageFilter->GetMedian(labelValue) << std::endl;
		//			std::cout << "mean: " << labelStatisticsImageFilter->GetMean(labelValue) << std::endl;
		//			std::cout << "sigma: " << labelStatisticsImageFilter->GetSigma(labelValue) << std::endl;
		//			std::cout << "variance: " << labelStatisticsImageFilter->GetVariance(labelValue) << std::endl;
		//			std::cout << "sum: " << labelStatisticsImageFilter->GetSum(labelValue) << std::endl;
		//			std::cout << "count: " << labelStatisticsImageFilter->GetCount(labelValue) << std::endl;
		//			// std::cout << "box: " << labelStatisticsImageFilter->GetBoundingBox( labelValue ) << std::endl; // can't output
		//			// a box
		//			std::cout << "region: " << labelStatisticsImageFilter->GetRegion(labelValue) << std::endl;
		//			std::cout << std::endl << std::endl;
		//		}
		//	}

		

			//rozrost 
		/*	ConnectedThreshold::Pointer connThres = ConnectedThreshold::New();
			connThres->SetInput(image3D);
			connThres->SetConnectivity(ConnectedThreshold::FullConnectivity);
			connThres->SetLower(500);
		connThres->SetUpper(200);
			ConnectedThreshold::IndexType seed1;
			seed1[0] = 129; seed1[1] = 136; seed1[2] = 48;
			connThres->SetSeed(seed1);
			connThres->Update();

			series3DWriter->SetInput(connThres->GetOutput());
			series3DWriter->SetFileName("..\\wyniki\\rozrost.vtk");
			series3DWriter->Update();*/
			///spróbuj Confidence Filter filtr medianowy 


			//otwarcie
			/*BallType::SizeType rad;
			rad[0] = 9;
			rad[1] = 5;
			int radius = 1;
			using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, Image3DType::ImageDimension>;
			StructuringElementType structuringElement;
			structuringElement.SetRadius(radius);
			structuringElement.CreateStructuringElement();

			using BinaryMorphologicalOpeningImageFilterType = itk::BinaryMorphologicalOpeningImageFilter <Image3DType, Image3DType, StructuringElementType>;
			BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter
				= BinaryMorphologicalOpeningImageFilterType::New();
			openingFilter->SetInput(thresholder->GetOutput());
			openingFilter->SetKernel(structuringElement);
			openingFilter->SetForegroundValue(1);
			openingFilter->Update();

			series3DWriter->SetInput(openingFilter->GetOutput());
			series3DWriter->SetFileName("..\\wyniki\\img3D_op_seria.vtk");
			series3DWriter->Update();

			CCImageFilter->SetInput(image3D);
			CCImageFilter->Update();
			std::cout << CCImageFilter->GetObjectCount() << std::endl;*/


			//Magnitude

			//using FilterGradientType = itk::GradientMagnitudeImageFilter<Image3DType, Image3DType >;
			//FilterGradientType::Pointer filterG = FilterGradientType::New();
			//filterG->SetInput(image3D);
			//filterG->Update();

			//series3DWriter->SetInput(filterG->GetOutput());
			//series3DWriter->SetFileName("..\\wyniki\\img3D_mag.vtk");
			//series3DWriter->Update();


			//thresholder->SetInput(filterG->GetOutput());
			//thresholder->SetInsideValue(1);
			//thresholder->SetOutsideValue(0);
			//thresholder->SetLowerThreshold(120);
			//thresholder->SetUpperThreshold(180);

			//series3DWriter->SetInput(thresholder->GetOutput());
			//series3DWriter->SetFileName("..\\wyniki\\img3D_bin_po_mag.vtk");
			//series3DWriter->Update();




			//series3DWriter->SetInput(image3D);

			//writer->SetInput(openingFilter->GetOutput());
			//writer->SetFileName("../wyniki/otwarcie_48.dcm");
			//writer->Update();


			//using FilterType = itk::OtsuMultipleThresholdsImageFilter<ImageType, ImageType>;
			//FilterType::Pointer thresholderOtsu = FilterType::New();
			//thresholderOtsu->SetInput(image);
			//thresholderOtsu->SetLabelOffset(0);
			//thresholderOtsu->SetNumberOfHistogramBins(100);
			//thresholderOtsu->SetNumberOfThresholds(5);

			//writer->SetInput(thresholderOtsu->GetOutput());
			//writer->SetFileName("../wyniki/intensywnosc_48.dcm");
			//writer->Update();
	}
catch(itk::ExceptionObject &ex){
ex.Print(std::cout);
}
catch(std::exception &ex){
std::cout << ex.what() << std::endl;
}
catch(...){
std::cout << "Unknown error!" << std::endl;
}
std::cout << "Hit [Enter]...";
std::cin.get();
return EXIT_SUCCESS; // albo EXIT_FAILURE
}
