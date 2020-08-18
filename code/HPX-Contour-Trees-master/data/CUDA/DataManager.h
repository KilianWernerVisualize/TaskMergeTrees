#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <limits>
#include <string>
#include <vtkActor.h>
#include <vtkBYUReader.h>
#include <vtkExtractEdges.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtksys/SystemTools.hxx>

#include <vtkObjectFactory.h>
#include <vtkNrrdReader.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataSetMapper.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkImageDataToPointSet.h>
#include <vtkImageResample.h>
#include <vtkMarchingCubes.h>
#include <vtkPiecewiseFunction.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSampleFunction.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkStructuredData.h>
#include <vtkStructuredGrid.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>

class Options {
public:
    bool visualize;
    bool trunkSkipping;
    std::string inputFile;
    int resampleX;
    int resampleY;
    int resampleZ;
};

#include "TreeConstructor.h"

namespace std {
template <>
class numeric_limits<ctValue> {
public:
    static ctValue min()
    {
        ctValue result;
        result.value = std::numeric_limits<double>::min();
        result.key = std::numeric_limits<std::uint64_t>::min();
        return result;
    }

    static ctValue max()
    {
        ctValue result;
        result.value = std::numeric_limits<double>::max();
        result.key = std::numeric_limits<std::uint64_t>::max();
        return result;
    }
};
}

class Datamanager {

public:
    Datamanager()
    {
        this->resampleX = 100;
        this->resampleY = 100;
        this->resampleZ = 100;
    }

    Datamanager(int resampleX, int resampleY, int resampleZ)
    {
        this->resampleX = resampleX;
        this->resampleY = resampleY;
        this->resampleZ = resampleZ;
    }

    void readFromFile(const std::string& fileName)
    {
        std::cout << "Input: " << fileName << std::endl;

        std::string extension = vtksys::SystemTools::GetFilenameExtension(fileName);
        if (extension == ".vtu") {
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(fileName.c_str());
            reader->Update();
            extractEdges = vtkSmartPointer<vtkExtractEdges>::New();
            extractEdges->SetInputConnection(reader->GetOutputPort());
            extractEdges->Update();
            mesh = extractEdges->GetOutput();
            mesh->ComputeBounds();
            double bounds[6];
            mesh->GetBounds(bounds);
            min_value = bounds[0];
            max_value = bounds[1];
            conversionFactor = (std::numeric_limits<std::uint32_t>().max() / (max_value - min_value));

            vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
            mapper->SetInputConnection(extractEdges->GetOutputPort());

            actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
        } else if (extension == ".vti") {
            vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
            reader->SetFileName(fileName.c_str());
            reader->Update();

            std::cout << "Range: " << reader->GetOutput()->GetScalarRange()[0] << " - " << reader->GetOutput()->GetScalarRange()[1] << std::endl;

            vtkSmartPointer<vtkImageResample>
                resampler
                = vtkSmartPointer<vtkImageResample>::New();

            resampler->SetInputConnection(reader->GetOutputPort());
            resampler->SetOutputExtentToDefault();
            resampler->SetOutputOriginToDefault();
            resampler->SetAxisMagnificationFactor(0, resampleX / 100.0);
            resampler->SetAxisMagnificationFactor(1, resampleY / 100.0);
            resampler->SetAxisMagnificationFactor(2, resampleZ / 100.0);
            resampler->SetInterpolationModeToLinear();
            resampler->Update();

            grid = resampler->GetOutput();
            double* tmp = grid->GetScalarRange();
            min_value = tmp[0];
            max_value = tmp[1];
            conversionFactor = (std::numeric_limits<std::uint32_t>().max() / (max_value - min_value));
            values = static_cast<double*>(grid->GetScalarPointer());
            fvalues = (float*)malloc(grid->GetNumberOfPoints()*sizeof(float));
            for (int i = 0; (i < grid->GetNumberOfPoints()); i++){
                fvalues[i] = values[i];
            }
            grid->GetDimensions(dimensions);

        } else if (extension == ".nrrd") {
            vtkObjectFactory::SetAllEnableFlags(false, "vtkNrrdReader", "vtkPNrrdReader");
            vtkSmartPointer<vtkNrrdReader> reader = vtkSmartPointer<vtkNrrdReader>::New();
            reader->SetFileName(fileName.c_str());
            reader->Update();

            std::cout << "Range: " << reader->GetOutput()->GetScalarRange()[0] << " - " << reader->GetOutput()->GetScalarRange()[1] << std::endl;

            vtkSmartPointer<vtkImageResample>
                resampler
                = vtkSmartPointer<vtkImageResample>::New();

            resampler->SetInputConnection(reader->GetOutputPort());
            resampler->SetOutputExtentToDefault();
            resampler->SetOutputOriginToDefault();
            resampler->SetAxisMagnificationFactor(0, resampleX / 100.0);
            resampler->SetAxisMagnificationFactor(1, resampleY / 100.0);
            resampler->SetAxisMagnificationFactor(2, resampleZ / 100.0);
            resampler->SetInterpolationModeToLinear();
            resampler->Update();

            grid = resampler->GetOutput();
            double* tmp = grid->GetScalarRange();
            min_value = tmp[0];
            max_value = tmp[1];
            conversionFactor = (std::numeric_limits<std::uint32_t>().max() / (max_value - min_value));
            fvalues = static_cast<float*>(grid->GetScalarPointer());
            grid->GetDimensions(dimensions);
            std::cout << "Dimensions: " << dimensions[0] << "/" << dimensions[1] << "/" << dimensions[2];
        } else {
            std::cout << "Error: Unknown file format" << std::endl;
        }

        // Precompute neighbors
        this->neighborsBuffer.reserve(getNumVertices() * 6);
        this->neighborsMap.reserve(getNumVertices());

        std::vector<int> neighbors;

        for (int i = 0; i < getNumVertices(); i++) {
            this->computeNeighbors(i, neighbors);

            NeighborsEntry e;
            e.numNeighbors = neighbors.size();
            e.offset = neighborsBuffer.size();

            this->neighborsMap.push_back(e);

            for (int n : neighbors)
                this->neighborsBuffer.push_back(n);
        }

        std::cout << "Vertices: " << this->getNumVertices() << std::endl;
    }

    const int* getNeighbors(int v, int* numNeighborsOut) const
    {
        *numNeighborsOut = this->neighborsMap[v].numNeighbors;
        return &this->neighborsBuffer[this->neighborsMap[v].offset];
    }

    ctValue getValue(int v) const
    {
        if (mesh != nullptr) {
            double p[3];
            mesh->GetPoint(v, p);
            ctValue result;
            result.value = p[0];
            result.key = v;
            return result;
        } else {
            ctValue result;
            result.value = values[v];
            result.key = v;
            return result;
        }
    }

    void getVertex(int v, double* p) const
    {
        if (mesh != nullptr)
            mesh->GetPoint(v, p);
        else
            grid->GetPoint(v, p);
    }

    int getNumVertices() const
    {
        if (mesh != nullptr)
            return (int) mesh->GetNumberOfPoints();
        else
            return (int) grid->GetNumberOfPoints();
    }

    void render(vtkSmartPointer<vtkRenderer> renderer)
    {
        if (mesh != nullptr) {
            actor->GetProperty()->SetOpacity(0.5);
            actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
            renderer->AddActor(actor);
        } else {
            vtkSmartPointer<vtkSmartVolumeMapper> mapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
            mapper->SetBlendModeToComposite();
            mapper->SetInputData(grid);

            vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
            volumeProperty->ShadeOff();
            volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);

            vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
            compositeOpacity->AddPoint(0.0, 0.0);
            compositeOpacity->AddPoint(80.0, 0.1);
            compositeOpacity->AddPoint(80.1, 0.0);
            compositeOpacity->AddPoint(255.0, 0.0);
            volumeProperty->SetScalarOpacity(compositeOpacity); // composite first.

            vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
            color->AddRGBPoint(0.0, 0.0, 0.0, 1.0);
            color->AddRGBPoint(40.0, 1.0, 0.0, 0.0);
            color->AddRGBPoint(255.0, 1.0, 1.0, 1.0);
            volumeProperty->SetColor(color);

            volume = vtkSmartPointer<vtkVolume>::New();
            volume->SetMapper(mapper);
            volume->SetProperty(volumeProperty);

            vtkSmartPointer<vtkMarchingCubes> contour = vtkSmartPointer<vtkMarchingCubes>::New();
            contour->SetInputData(grid);
            contour->SetValue(0, 50);

            vtkSmartPointer<vtkPolyDataMapper> contourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            contourMapper->SetInputConnection(contour->GetOutputPort());
            contourMapper->ScalarVisibilityOff();

            actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(contourMapper);
            renderer->AddActor(actor);
            //renderer->AddViewProp(volume);
            renderer->ResetCamera();
        }
    }

    void printTimes()
    {
    }

    vtkSmartPointer<vtkActor> getRenderActor()
    {
        return actor;
    }

    double* rawValues(){
        return values;
    }

    std::vector<int>* rawNeighbors(){
        return &neighborsBuffer;
    }

    std::vector<NeighborsEntry>* neighborMap(){
        return &neighborsMap;
    }

    vtkPolyData* mesh = nullptr;
private:
    vtkSmartPointer<vtkExtractEdges> extractEdges;

    //    std::vector<std::vector<int>> allNeighbors;

    std::vector<NeighborsEntry> neighborsMap;
    std::vector<int> neighborsBuffer;

    vtkSmartPointer<vtkImageData> grid;
    int resampleX;
    int resampleY;
    int resampleZ;
    //vtkDataArray* values;
public:
    float* fvalues;
    double* values;
    double min_value;
    double max_value;
    double conversionFactor;

    mutable int dimensions[3];
private:
    vtkSmartPointer<vtkActor> actor;
    vtkSmartPointer<vtkVolume> volume;

    double factor = 1.0;

    vtkIdType coordToId(int ijk[]) const
    {
        return vtkStructuredData::ComputePointId(dimensions, ijk);
    }

    void idToCoord(vtkIdType id, int ijk[]) const
    {
        //vtkStructuredData::ComputePointStructuredCoords(id, dimensions, ijk);
        ijk[0] = id % dimensions[0];
        ijk[1] = (int)( id / dimensions[0] ) % dimensions[1];
        ijk[2] = (int)( id / ( dimensions[0]*dimensions[1] ));
    }

    vtkSmartPointer<vtkIdList> GetConnectedVertices(int id) const
    {
        vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();

        //get all cells that vertex 'id' is a part of
        vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
        mesh->GetPointCells(id, cellIdList);

        for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++) {
            vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
            mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

            if (pointIdList->GetId(0) != id)
                connectedVertices->InsertNextId(pointIdList->GetId(0));
            else
                connectedVertices->InsertNextId(pointIdList->GetId(1));
        }

        return connectedVertices;
    }

    void computeNeighbors(int v, std::vector<int>& nbrsBuf) const
    {
        nbrsBuf.clear();

        if (mesh != nullptr) {

            vtkSmartPointer<vtkIdList> connectedVertices = GetConnectedVertices(v);

            for (vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++) {
                nbrsBuf.push_back(connectedVertices->GetId(i));
            }

        } else {
            int coords[3];
            int neighbor[3];
            idToCoord(v, coords);
            nbrsBuf.reserve(6);
            if (coords[0] < dimensions[0] - 1) {
                neighbor[0] = coords[0] + 1;
                neighbor[1] = coords[1];
                neighbor[2] = coords[2];
                nbrsBuf.push_back(coordToId(neighbor));
            }
            if (coords[0] > 0) {
                neighbor[0] = coords[0] - 1;
                neighbor[1] = coords[1];
                neighbor[2] = coords[2];
                nbrsBuf.push_back(coordToId(neighbor));
            }
            if (coords[1] < dimensions[1] - 1) {
                neighbor[1] = coords[1] + 1;
                neighbor[0] = coords[0];
                neighbor[2] = coords[2];
                nbrsBuf.push_back(coordToId(neighbor));
            }
            if (coords[1] > 0) {
                neighbor[1] = coords[1] - 1;
                neighbor[0] = coords[0];
                neighbor[2] = coords[2];
                nbrsBuf.push_back(coordToId(neighbor));
            }
            if (coords[2] < dimensions[2] - 1) {
                neighbor[2] = coords[2] + 1;
                neighbor[1] = coords[1];
                neighbor[0] = coords[0];
                nbrsBuf.push_back(coordToId(neighbor));
            }
            if (coords[2] > 0) {
                neighbor[2] = coords[2] - 1;
                neighbor[1] = coords[1];
                neighbor[0] = coords[0];
                nbrsBuf.push_back(coordToId(neighbor));
            }
        }
    }
};

#endif
