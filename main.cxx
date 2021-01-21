#include <iostream> // std::cout, std::cin
#include <string> // std::string
#include <itkImage.h> // itk::Image
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageSeriesReader.h>
#include <itkImageSeriesWriter.h>
#include <itkNumericSeriesFileNames.h>
#include <itkGDCMSeriesFileNames.h> 
#include <itkGDCMImageIO.h>
#include<itkIntensityWindowingImageFilter.h>
#include<itkThresholdImageFilter.h>
#include<itkBinaryThresholdImageFilter.h>
#include<itkBinaryBallStructuringElement.h>
#include<itkErodeObjectMorphologyImageFilter.h>
#include<itkBinaryErodeImageFilter.h>
#include<itkBinaryDilateImageFilter.h>
#include<itkConnectedComponentImageFilter.h>
#include<itkConnectedThresholdImageFilter.h>
#include<itkMultiplyImageFilter.h>
#include<itkBinaryMorphologicalOpeningImageFilter.h>
#include<itkConfidenceConnectedImageFilter.h>
#include <itkInvertIntensityImageFilter.h>
#include<itkRelabelComponentImageFilter.h>
#include<itkFlipImageFilter.h>
#include<itkAddImageFilter.h>
#include "itkLabelGeometryImageFilter.h"





using PixelType = signed short;
using ImageType = itk::Image<PixelType, 2>;
using BallType = itk::BinaryBallStructuringElement<short, 3>;


using Image3DType = itk::Image<PixelType, 3>;
using ImageIOType = itk::GDCMImageIO;
//###################### Odczytywanie i zapis plików #################
using SeriesReaderType = itk::ImageSeriesReader<Image3DType>;
using Series2DWriterType = itk::ImageSeriesWriter<Image3DType, ImageType>;
using Series3DWriterType = itk::ImageSeriesWriter<Image3DType, Image3DType>;
using Writer3Dtype = itk::ImageFileWriter<Image3DType>;
using ReaderType = itk::ImageFileReader<ImageType>; 
using Reader3DType = itk::ImageFileReader<Image3DType>;
using WriterType = itk::ImageFileWriter<ImageType>;
//#####################################################################
using ConnectedComponent = itk::ConnectedComponentImageFilter<Image3DType, Image3DType, Image3DType>;
using NumSeriesFileNames = itk::NumericSeriesFileNames;
using ConnectedThreshold = itk::ConnectedThresholdImageFilter<Image3DType, Image3DType>;
using FilterType = itk::BinaryThresholdImageFilter<Image3DType, Image3DType>;
using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, Image3DType::ImageDimension>;
using FilterDilateType = itk::BinaryDilateImageFilter<Image3DType, Image3DType, StructuringElementType>;
using FilterErodeType = itk::BinaryErodeImageFilter<Image3DType, Image3DType, StructuringElementType>;
using MultiplyType = itk::MultiplyImageFilter< Image3DType, Image3DType, Image3DType>;
using OpeningFilterType = itk::BinaryMorphologicalOpeningImageFilter<Image3DType, Image3DType, StructuringElementType>;
using ConnCompFilterType = itk::ConnectedComponentImageFilter<Image3DType, Image3DType>;
using LabelGeometryImageFilterType = itk::LabelGeometryImageFilter<Image3DType>;
using ConnectedFilterType = itk::ConfidenceConnectedImageFilter<Image3DType, Image3DType>;
using InvertIntensityImageFilterType = itk::InvertIntensityImageFilter<Image3DType>;
using WindowingImageFilter = itk::IntensityWindowingImageFilter<Image3DType>;
using RelabelComponentFilterType = itk::RelabelComponentImageFilter<Image3DType, Image3DType>;
using TransformType = itk::AffineTransform<double, 3>;
using AddImageFilterType = itk::AddImageFilter<Image3DType, Image3DType>;

int main() // glowna funkcja programu
{
	try {
		ReaderType::Pointer reader = ReaderType::New();
		Reader3DType::Pointer reader3D = Reader3DType::New();
		Reader3DType::Pointer reader3Db = Reader3DType::New();
		WriterType::Pointer writer = WriterType::New();
		ImageType::Pointer image = ImageType::New();

		SeriesReaderType::Pointer seriesReader = SeriesReaderType::New();
		Series2DWriterType::Pointer series2DWriter = Series2DWriterType::New();
		itk::GDCMSeriesFileNames::Pointer namesGen = itk::GDCMSeriesFileNames::New();
		NumSeriesFileNames::Pointer numSeriesFileNames = NumSeriesFileNames::New();
		Image3DType::Pointer image3D = Image3DType::New();
		Image3DType::Pointer image3D_bin = Image3DType::New();
		Image3DType::Pointer image3D_mnozenie = Image3DType::New();
		Image3DType::Pointer image3D_mnozenie_left = Image3DType::New();

		ImageIOType::Pointer dicomIO = ImageIOType::New();
		Writer3Dtype::Pointer writer3D = Writer3Dtype::New();
		Series3DWriterType::Pointer series3DWriter = Series3DWriterType::New();
		ConnectedComponent::Pointer CCImageFilter = ConnectedComponent::New();
		itk::GDCMImageIO::Pointer gdcmImageIO = itk::GDCMImageIO::New();

		namesGen->SetDirectory("../dane_glowa/glowa/mozg zb/Head_Neck_Standard - 153275/t2_tse_tra_5");
		itk::GDCMSeriesFileNames::SeriesUIDContainerType seriesUIDs = namesGen->GetSeriesUIDs();
		itk::GDCMSeriesFileNames::FileNamesContainerType fileNames = namesGen->GetFileNames(seriesUIDs[0]);
		for (size_t i = 0; i < seriesUIDs.size(); i++) {
			std::cout << "UID serii:"+ seriesUIDs[i] << std::endl;
		}
		//######################## Czytanie serii obrazowej ##########################
		std::cout << "Czytanie i zapisywanie serii obrazowej..." << std::endl;

		seriesReader->SetFileNames(namesGen->GetFileNames(seriesUIDs[0]));
		seriesReader->SetImageIO(gdcmImageIO);
		image3D = seriesReader->GetOutput();
		seriesReader->Update();


		//############### Zapisywanie serii obrazowej jako pliki dicom #############################

		//numSeriesFileNames->SetSeriesFormat("..\\wyniki\\IMG\%05d.dcm");
		//numSeriesFileNames->SetStartIndex(1);
		//numSeriesFileNames->SetEndIndex(image3D->GetLargestPossibleRegion().GetSize()[2]);

		//series2DWriter->SetFileNames(numSeriesFileNames->GetFileNames());
	    //series2DWriter->SetImageIO(gdcmImageIO);
		////series2DWriter->SetMetaDataDictionaryArray(seriesReader->GetMetaDataDictionaryArray());
		//series2DWriter->SetInput(image3D);
		//series2DWriter->Update();
		
		//########## Zapisywanie obrazu jako jeden plik 3D z rozszerzeniem .vtk ####################
		series3DWriter->SetInput(image3D);
		series3DWriter->SetFileName("..\\wyniki\\1_img3D.vtk");
		series3DWriter->Update();


		//########### Binaryzacja ######################
		std::cout << "Binaryzacja..." << std::endl;

		FilterType::Pointer thresholder = FilterType::New();
		thresholder->SetInput(image3D);
		thresholder->SetInsideValue(1);
		thresholder->SetOutsideValue(0);
		thresholder->SetLowerThreshold(300);
		series3DWriter->SetInput(thresholder->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\2_img3D_bin_calosc.vtk");
		series3DWriter->Update();
		image3D_bin = thresholder->GetOutput();
		image3D_bin->CopyInformation(image3D);


		//########### Dylatacja ########################
		std::cout << "Dylatacja..." << std::endl;

		BallType::SizeType rad;
		rad[0] = 8;
		rad[1] = 8;
		rad[2] = 8;
		StructuringElementType structuringElement;
		structuringElement.SetRadius(rad);
		structuringElement.CreateStructuringElement();

		FilterDilateType::Pointer dilateFilter = FilterDilateType::New();
		dilateFilter->SetInput(image3D_bin);
		dilateFilter->SetKernel(structuringElement);
		dilateFilter->SetBackgroundValue(0);
		dilateFilter->SetForegroundValue(1);
		dilateFilter->Update();

		series3DWriter->SetInput(dilateFilter->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\3_img3D_dilate.vtk");
		series3DWriter->Update();

		//########### Erozja ############################# 
		std::cout << "Erozja..." << std::endl;

		rad[0] = 38;
		rad[1] = 38;
		rad[2] = 38;
		structuringElement.SetRadius(rad);
		structuringElement.CreateStructuringElement();

		FilterErodeType::Pointer erodeFilter = FilterErodeType::New();
		erodeFilter->SetInput(dilateFilter->GetOutput());
		erodeFilter->SetKernel(structuringElement);
		erodeFilter->SetBackgroundValue(0);
		erodeFilter->SetForegroundValue(1);
		erodeFilter->Update();

		series3DWriter->SetInput(erodeFilter->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\4_img3D_erozja.vtk");
		series3DWriter->Update();

		//############# Mno¿enie maski z obrazem oryginalnym ##############
		std::cout << "Mno¿enie maski z obrazem oryginalnym..." << std::endl;

		MultiplyType::Pointer multiply = MultiplyType::New();
		multiply->SetInput1(image3D);
		multiply->SetInput2(erodeFilter->GetOutput());
		image3D_mnozenie = multiply->GetOutput();
		image3D_mnozenie->CopyInformation(image3D);
		multiply->Update();
		series3DWriter->SetInput(image3D_mnozenie);
		series3DWriter->SetFileName("..\\wyniki\\5_img3D_mnozenie.vtk");
		series3DWriter->Update();

		

		//############ Binaryzacja ##########################################
		std::cout << "Binaryzacja..." << std::endl;

		thresholder->SetInput(multiply->GetOutput());
		thresholder->SetInsideValue(1);
		thresholder->SetOutsideValue(0);
		thresholder->SetLowerThreshold(600);
		thresholder->SetUpperThreshold(800);
		series3DWriter->SetInput(thresholder->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\6_img3D_bin_komory.vtk");
		series3DWriter->Update();

		//############# Otwarcie #####################
		std::cout << "Otwarcie..." << std::endl;

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
		series3DWriter->SetFileName("..\\wyniki\\7_img3D_otwarcie.vtk");
		series3DWriter->Update();

		//################# Etykietowanie ################
		std::cout << "Etykietowanie..." << std::endl;

		ConnCompFilterType::Pointer connComp1 = ConnCompFilterType::New();
		connComp1->SetInput(openFilter->GetOutput());
		connComp1->SetBackgroundValue(0);
		connComp1->Update();
		series3DWriter->SetInput(connComp1->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\8_img3D_etykiety.vtk");
		series3DWriter->Update();

		//############### Pobieranie centroidu ostatniej etykiety (odpowiadaj¹cej czêœci komory prawa) ######################
		std::cout << "Pobieranie centroidu ostatniej etykiety komory..." << std::endl;

		LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
		labelGeometryImageFilter->SetInput(connComp1->GetOutput());
		labelGeometryImageFilter->CalculatePixelIndicesOn();
		labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();
		labelGeometryImageFilter->Update();
		LabelGeometryImageFilterType::LabelsType labelsAll= labelGeometryImageFilter->GetLabels();
		LabelGeometryImageFilterType::LabelPointType startPoint = labelGeometryImageFilter->GetCentroid(labelsAll.size() - 1);
		std::cout << "Wspolrzedne centroidu:" << std::endl;
		std::cout << startPoint[0] << std::endl;
		std::cout << startPoint[1] << std::endl;
		std::cout << startPoint[2] << std::endl;


		//################### Rozrost obszaru od punktu wyznaczonego w poprzednim kroku ###############################
		std::cout << "Rozrost obszaru..." << std::endl;

		ConnectedFilterType::Pointer confidenceConnected = ConnectedFilterType::New();
		confidenceConnected->SetInput(image3D); 
		confidenceConnected->SetMultiplier(1.3); 
		confidenceConnected->SetNumberOfIterations(1);
		confidenceConnected->SetReplaceValue(99);
		confidenceConnected->SetInitialNeighborhoodRadius(2);
		ConnectedFilterType::IndexType index;
		index[0] = (int)(startPoint[0]) + 2; 
		index[1] = (int)(startPoint[1])+2; 
		index[2] = (int)(startPoint[2])+2;

		confidenceConnected->AddSeed(index);
		confidenceConnected->Update();
		series3DWriter->SetInput(confidenceConnected->GetOutput());
		series3DWriter->SetFileName("..\\wyniki\\9_img3D_rozrost-thresh.vtk");
		series3DWriter->Update();
	//##################### Odwracanie intensywnoœci ###################
		std::cout << "Odwracanie intensywnosci..." << std::endl;

	InvertIntensityImageFilterType::Pointer invertIntensityFilter = InvertIntensityImageFilterType::New();
	invertIntensityFilter->SetInput(confidenceConnected->GetOutput());
	invertIntensityFilter->SetMaximum(1);
	invertIntensityFilter->Update();
	series3DWriter->SetInput(invertIntensityFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9a_img3D_invert.vtk");
	series3DWriter->Update();

	//################ Okno intensywnoœci #################################
	std::cout << "Okno intensywnosci..." << std::endl;

	WindowingImageFilter::Pointer windowingImageFilter = WindowingImageFilter::New();
	windowingImageFilter->SetInput(invertIntensityFilter->GetOutput());
	windowingImageFilter->SetWindowMaximum(1);
	windowingImageFilter->SetWindowMinimum(-1);
	windowingImageFilter->SetOutputMinimum(0);
	windowingImageFilter->SetOutputMaximum(1);
	windowingImageFilter->Update();
	series3DWriter->SetInput(windowingImageFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9b_img3D_window.vtk");
	series3DWriter->Update();
	
	Image3DType::Pointer windowImage = Image3DType::New();
	Image3DType::Pointer erode = Image3DType::New();


	//################## Mno¿enie obrazów ######################
	std::cout << "Mnozenie obrazow..." << std::endl;

	reader3D->SetFileName("..\\wyniki\\9b_img3D_window.vtk");
	reader3D->Update();
	windowImage = reader3D->GetOutput();

	reader3Db->SetFileName("..\\wyniki\\4_img3D_erozja.vtk");
	reader3Db->Update();
	erode = reader3Db->GetOutput();
	//erode->CopyInformation(windowImage);

	multiply->SetInput1(windowImage);
	multiply->SetInput2(erode);
	multiply->Update();
	series3DWriter->SetInput(multiply->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9c_img3D_mnozenie_komory.vtk");
	series3DWriter->Update();


	//################################ Etykietowanie ##########################
	std::cout << "Etykietowanie..." << std::endl;

	connComp1->SetInput(multiply->GetOutput());
	connComp1->SetBackgroundValue(0);
	connComp1->Update();

	series3DWriter->SetInput(connComp1->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9d_img3D_label.vtk");
	series3DWriter->Update();

	//################### Reetykietowanie, usuniêcie obiektów mniejszych ni¿ 100 pikseli ##############################
	std::cout << "Reetykietowanie..." << std::endl;

	RelabelComponentFilterType::Pointer relabel = RelabelComponentFilterType::New();
	relabel->SetInput(connComp1->GetOutput());
	relabel->SetMinimumObjectSize(100);
	
	relabel->Update();

	series3DWriter->SetInput(relabel->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9e_img3D_relabel.vtk");
	series3DWriter->Update();

	//######################### Wybór etykiety -- wysegmentowana prawa komora #################################
	std::cout << "Wybór etykiety..." << std::endl;

	FilterType::Pointer thresholder2 = FilterType::New();
	thresholder2->SetInput(relabel->GetOutput());
	thresholder2->SetLowerThreshold(2);
	thresholder2->SetUpperThreshold(2);
	thresholder2->Update();

	series3DWriter->SetInput(thresholder2->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9f_img3D_binaryzacja1komora.vtk");
	series3DWriter->Update();
	
	//##################### Odbicie maski prawej komory na lew¹ pó³kule ###################
	std::cout << "Odbicie maski..." << std::endl;

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
	series3DWriter->SetFileName("..\\wyniki\\9g_img3D_flip.vtk");
	series3DWriter->Update();
	Image3DType::Pointer flipImage = Image3DType::New();

	flipImage = flipFilter->GetOutput();
	flipImage->CopyInformation(thresholder2->GetOutput());
	
	//#################### Dylatacja obszaru poszukiwania maski lewej ######################
	std::cout << "Dylatacja obszaru poszukiwania maski lewej..." << std::endl;

	FilterDilateType::Pointer dilateMask = FilterDilateType::New();

	windowingImageFilter->SetInput(flipImage);
	windowingImageFilter->SetWindowMaximum(1);
	windowingImageFilter->SetWindowMinimum(0);
	windowingImageFilter->SetOutputMinimum(0);
	windowingImageFilter->SetOutputMaximum(1);
	windowingImageFilter->Update();

	StructuringElementType dilateElement;
	rad[0] = 8;
	rad[1] = 3;
	rad[2] = 1;
	dilateElement.SetRadius(rad);
	dilateElement.CreateStructuringElement();

	dilateFilter->SetInput(windowingImageFilter->GetOutput());
	dilateFilter->SetKernel(dilateElement);
	dilateFilter->SetBackgroundValue(0);
	dilateFilter->SetForegroundValue(1);
	dilateFilter->Update();

	series3DWriter->SetInput(dilateFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9h_img3D_dilate.vtk");
	series3DWriter->Update();


	//#####################	Mno¿enie prawej przesuniêtej maski z lew¹ komor¹  ######################
	
	MultiplyType::Pointer multiply_labelled = MultiplyType::New();
	multiply_labelled->SetInput1(relabel->GetOutput());
	image3D_mnozenie_left = dilateFilter->GetOutput();
	image3D_mnozenie_left->CopyInformation(relabel->GetOutput());
	multiply_labelled->SetInput2(image3D_mnozenie_left);

	multiply_labelled->Update();
	series3DWriter->SetInput(multiply_labelled->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9i_img3D_left_v.vtk");
	series3DWriter->Update();
//############### Po³¹cznie maski lewej i prawej komory ########################
	std::cout << "Polaczenie masek..." << std::endl;

	AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

	Image3DType::Pointer left_v = Image3DType::New();
	left_v = multiply_labelled->GetOutput();
	left_v->CopyInformation(thresholder2->GetOutput());


	addFilter->SetInput1(left_v);
	addFilter->SetInput2(thresholder2->GetOutput());


	series3DWriter->SetInput(addFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9j_img3D_dodawanie.vtk");
	series3DWriter->Update();

	//################# Okno intensywnoœci ######################33
	std::cout << "Okno intensywnosci..." << std::endl;

	WindowingImageFilter::Pointer windowingFilter = WindowingImageFilter::New();
	windowingFilter->SetInput(addFilter->GetOutput());
	windowingFilter->SetWindowMaximum(1);
	windowingFilter->SetWindowMinimum(0);
	windowingFilter->SetOutputMinimum(0);
	windowingFilter->SetOutputMaximum(1000);
	windowingFilter->Update();
	series3DWriter->SetInput(windowingFilter->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\9k_img3D_dodawanie_window.vtk");
	series3DWriter->Update();

	//#################### Na³o¿enie masek komór na obraz oryginalny ####################
	std::cout << "Nalozenie masek na komor na obrazy oryginalne..." << std::endl;

	Image3DType::Pointer maskImage = Image3DType::New();
	maskImage = windowingFilter->GetOutput();
	maskImage->CopyInformation(image3D);

	AddImageFilterType::Pointer addFilter2 = AddImageFilterType::New();
	addFilter2->SetInput1(maskImage);
	addFilter2->SetInput2(image3D);


	series3DWriter->SetInput(addFilter2->GetOutput());
	series3DWriter->SetFileName("..\\wyniki\\Wynik_koncowy.vtk");
	series3DWriter->Update();





	std::cout << "Koniec algorytmu." << std::endl;



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
