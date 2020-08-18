#ifndef VISUALIZER_H
#define VISUALIZER_H

#include <hpx/hpx.hpp>

#include "Arc.h"
#include "DataManager.h"

#include <string>
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkExtractEdges.h>
#include <vtkExtractSelection.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkObjectFactory.h>
#include <vtkPlaneSource.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPropPicker.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>
#include <vtkTubeFilter.h>
#include <vtkVersion.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkPointSource.h>

//*****************************************   VTK Visualization *********************************************//
class MouseInteractorStyle2 : public vtkInteractorStyleTrackballCamera {
public:
    static MouseInteractorStyle2* New();
    vtkTypeMacro(MouseInteractorStyle2, vtkInteractorStyleTrackballCamera);

    virtual void OnLeftButtonDown()
    {
        int* clickPos = this->GetInteractor()->GetEventPosition();

        // Pick from this location.
        vtkSmartPointer<vtkPropPicker> picker = vtkSmartPointer<vtkPropPicker>::New();
        picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

        double* pos = picker->GetPickPosition();
        std::cout << "Pick position (world coordinates) is: "
                  << pos[0] << " " << pos[1]
                  << " " << pos[2] << std::endl;
        if (pos[0] == 0 && pos[1] == 0 && pos[2] == 0) {
            vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
            return;
        }
        std::cout << "Picked actor: " << picker->GetActor() << std::endl;
        if (picker->GetActor() == nullptr) {
            vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
            return;
        }
        if (picker->GetActor()->GetProperty()->GetOpacity() != 1.0) {
            vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
            return;
        }

        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

private:
    vtkSmartPointer<vtkActor> last = vtkSmartPointer<vtkActor>::New();
};

//*****************************************   VTK Visualization *********************************************//

class Visualizer {
    Datamanager* d;

public:
    std::vector<vtkSmartPointer<vtkActor>> minima_actors;
    std::vector<vtkSmartPointer<vtkActor>> branch_actors;
    std::vector<vtkSmartPointer<vtkActor>> boundary_actors;

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    hpx::lcos::local::mutex lock;

    void init(Datamanager& data)
    {
        d = &data;
    }

    void show()
    {

        vtkSmartPointer<vtkPolyData> pointsPolydata =
           vtkSmartPointer<vtkPolyData>::New();

         pointsPolydata->SetPoints(points);

         vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
           vtkSmartPointer<vtkVertexGlyphFilter>::New();

         vertexFilter->SetInputData(pointsPolydata);
         vertexFilter->Update();

         vtkSmartPointer<vtkPolyData> polydata =
           vtkSmartPointer<vtkPolyData>::New();
         polydata->ShallowCopy(vertexFilter->GetOutput());

         vtkSmartPointer<vtkPolyDataMapper> mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
          mapper->SetInputData(polydata);

          vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
          actor->SetMapper(mapper);
          actor->GetProperty()->SetPointSize(5);

        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);

        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        vtkSmartPointer<MouseInteractorStyle2> style = vtkSmartPointer<MouseInteractorStyle2>::New();
        style->SetDefaultRenderer(renderer);

        renderWindowInteractor->SetInteractorStyle(style);

        for (int i = 0; i < minima_actors.size(); i++)
            renderer->AddActor(minima_actors.at(i));
        for (int i = 0; i < branch_actors.size(); i++)
            renderer->AddActor(branch_actors.at(i));
        for (int i = 0; i < boundary_actors.size(); i++)
            renderer->AddActor(boundary_actors.at(i));

        renderer->AddActor(actor);
        renderer->SetBackground(.3, .6, .3);

        //d->render(renderer);

        renderWindow->Render();
        renderWindowInteractor->Initialize();
        renderWindowInteractor->Start();
    }

    void visualizePair(vtkIdType one, vtkIdType two)
    {
        return;
        double p[3];
        d->getVertex(one, p);

        vtkSmartPointer<vtkSphereSource> sphereSourceext = vtkSmartPointer<vtkSphereSource>::New();
        sphereSourceext->SetCenter(p);
        sphereSourceext->SetRadius(4);

        vtkSmartPointer<vtkPolyDataMapper> mapperext = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperext->SetInputConnection(sphereSourceext->GetOutputPort());

        vtkSmartPointer<vtkActor> actorext = vtkSmartPointer<vtkActor>::New();
        actorext->SetMapper(mapperext);
        actorext->GetProperty()->SetColor(1.0, 0.2, 0.4);
        actorext->GetProperty()->SetOpacity(0.4);

        minima_actors.push_back(actorext);

        vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
        lineSource->SetPoint1(p);

        double p2[3];
        d->getVertex(two, p2);

        p2[0] -= p[0];
        p2[1] -= p[1];
        p2[2] -= p[2];

        p2[0] /= 2;
        p2[1] /= 2;
        p2[2] /= 2;

        p2[0] += p[0];
        p2[1] += p[1];
        p2[2] += p[2];

        lineSource->SetPoint2(p2);

        vtkSmartPointer<vtkPolyDataMapper> lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        lineMapper->SetInputConnection(lineSource->GetOutputPort());
        vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
        lineActor->SetMapper(lineMapper);

        vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
        tubeFilter->SetInputConnection(lineSource->GetOutputPort());
        tubeFilter->SetRadius(0.05); //default is .5
        tubeFilter->SetNumberOfSides(50);
        tubeFilter->Update();

        // Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        tubeMapper->SetInputConnection(tubeFilter->GetOutputPort());
        vtkSmartPointer<vtkActor> tubeActor = vtkSmartPointer<vtkActor>::New();
        tubeActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
        if (one < 0 || two < 0)
            tubeActor->GetProperty()->SetColor(0.0, 1.0, 0.0);
        tubeActor->SetMapper(tubeMapper);

        branch_actors.push_back(tubeActor);

        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(p);
        sphereSource->SetRadius(4);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphereSource->GetOutputPort());

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1.0, 0.2, 0.4);
        actor->GetProperty()->SetOpacity(0.4);

        minima_actors.push_back(actor);
    }

    void visualizeBranch(Arc* newBranch, bool trunk)
    {

        double p[3];
        d->getVertex(newBranch->extremum, p);

        points->InsertNextPoint(p);

        /*
        vtkSmartPointer<vtkPointSource> pointSource = vtkSmartPointer<vtkPointSource>::New();
        pointSource->SetCenter(p);
        pointSource->SetRadius(1);
        pointSource->set

        vtkSmartPointer<vtkSphereSource> sphereSourceext = vtkSmartPointer<vtkSphereSource>::New();
        sphereSourceext->SetCenter(p);
        sphereSourceext->SetRadius(.04);

        vtkSmartPointer<vtkPolyDataMapper> mapperext = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperext->SetInputConnection(sphereSourceext->GetOutputPort());

        vtkSmartPointer<vtkActor> actorext = vtkSmartPointer<vtkActor>::New();
        actorext->SetMapper(mapperext);
        actorext->GetProperty()->SetColor(1.0, 0.2, 0.4);
        actorext->GetProperty()->SetOpacity(4);

        minima_actors.push_back(actorext);
        */

        vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
        lineSource->SetPoint1(p);
        d->getVertex(newBranch->saddle, p);
        lineSource->SetPoint2(p);
        lineSource->Update();

        vtkSmartPointer<vtkPolyDataMapper> lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        lineMapper->SetInputConnection(lineSource->GetOutputPort());
        vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
        lineActor->SetMapper(lineMapper);
        lineActor->GetProperty()->SetLineWidth(1);

        if (trunk)
            lineActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        else
            lineActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
        if (newBranch->extremum < 0 || newBranch->saddle < 0)
            lineActor->GetProperty()->SetColor(0.0, 1.0, 0.0);

        branch_actors.push_back(lineActor);

        points->InsertNextPoint(p);

        /*
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(p);
        sphereSource->SetRadius(1);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphereSource->GetOutputPort());

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1.0, 0.2, 0.4);
        actor->GetProperty()->SetOpacity(0.4);

        minima_actors.push_back(actor);
        */
    }

    void createBoundaryActor(int b, float r, float g, float blue, float alpha)
    {
        return;
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();

        double p[3];
        d->getVertex(b, p);

        sphereSource->SetCenter(p);
        sphereSource->SetRadius(0.2);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphereSource->GetOutputPort());

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(r, g, blue);
        actor->GetProperty()->SetOpacity(alpha);

        boundary_actors.push_back(actor);
    }

    void createMinimumActor(int v)
    {
        return;
        double p[3];
        d->getVertex(v, p);
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(p);
        sphereSource->SetRadius(0.4);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphereSource->GetOutputPort());

        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        actor->GetProperty()->SetOpacity(0.5);
        //if (swept[v] < 0)
        //    actor->GetProperty()->SetColor(0.0,1.0,0.0);

        int numNeighbors;
        d->getNeighbors(v, &numNeighbors);
        if (numNeighbors <= 0)
            actor->GetProperty()->SetColor(0.0, 1.0, 1.0);

        minima_actors.push_back(actor);
    }
};

#endif
