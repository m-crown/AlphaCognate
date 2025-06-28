import { useParams } from "react-router-dom";
import { useEffect, useState, useRef } from "react";
import TransplantTable from "../components/TransplantTable";
import { Grid, Loader, ScrollArea } from "@mantine/core";
import { StructureInfo } from "../components/StructureInfo";
import MolStarWrapper from "../components/MolViewSpec";
// import MolstarViewer from "../components/MolStarViewer";
// import CognateLigandsTable from "../components/CognateLigandsTable";

import pako from "pako";

export type CognateLigand = {
  id: number;
  name: string;
  smiles: string;
  xref: string;
  similarity: number;
};

export type RawTransplant = {
  transplant_id: number;
  structure_name: string;
  ligand_id: number;
  ligand_chain: string;
  ligand_residues: string;
  ec_list: string;
  tcs: number;
  foldseek_rmsd: number;
  global_rmsd: number;
  local_rmsd: number;
  hetcode: string;
  name: string;
  smiles: string;
};

export type TransplantRecord = {
  transplant: RawTransplant;
  cognate_ligands: CognateLigand[];
};

export type TransplantApiResponse = {
  data: TransplantRecord[];
  structure_data: Structure;
};

//this is used in multiple places - abstract away
export type Structure = {
  name: string;
  url: string;
  transplanted: boolean;
};

function StructurePage() {
  const { id } = useParams(); // Extract the ID from the URL
  const [transplants, setTransplants] = useState<TransplantRecord[]>([]);
  const [structure, setStructure] = useState<Structure | null>(null);
  const [loading, setLoading] = useState(true);
  //make a url for the structure viewer
  const [structure_url, setStructureUrl] = useState<string | null>(null);
  const [fileData, setFileData] = useState<string | null>(null);
  const viewerInstanceRef = useRef<any>(null);
  const [best, setBest] = useState(true);
  const viewerContainerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);

  useEffect(() => {

    let objectUrl: string | null = null;

    if (id) {
      const fetchTransplantData = async () => {
        const url = new URL(
          '/transplants',
          process.env.NODE_ENV === 'production'
            ? 'http://localhost:8000/' //replace with prod url when it exists
            : 'http://localhost:8000/',
        );
        url.searchParams.append('structure_name', id);
        if (best) url.searchParams.append('best', 'true');
        try {
          const response = await fetch(url.href);
          const json = (await response.json()) as TransplantApiResponse;
          if (!response.ok) {
            throw new Error('Failed to fetch data');
          }
          setTransplants(json.data); // Update state with fetched data
          setStructure(json.structure_data);
          setLoading(false); // Set loading to false after data is fetched
          return json.structure_data.name;
        } catch (error) {
          console.error(error);
          setLoading(false); // Set loading to false even if an error occurs
        }
      };

      const fetchStructureData = async (name: string) => {
        try {
          // download .gz file and decompress
          let structure_url = new URL(
            name,
            process.env.NODE_ENV === 'production'
              ? 'http://localhost:8000/structure/' //replace with prod url when it exists
              : 'http://localhost:8000/structure/',)
          // structure_url.searchParams.append('chains', "A");

          const response = await fetch(structure_url.href);
          // const response = await fetch(url);
          const compressedData = await response.arrayBuffer();
          
          const decompressedData = pako.ungzip(new Uint8Array(compressedData), { to: "string" });
          const blob = new Blob([decompressedData], { type: "text/plain" });

          // Create an object URL that points to it
          objectUrl = URL.createObjectURL(blob);
          console.log(objectUrl);
          // objectUrl = URL.createObjectURL(new Blob([decompressedData], { type: "text/plain" }));
          setStructureUrl(objectUrl);
        } catch (err) {
          console.error("Error fetching or decompressing file:", err);
        }
      };
      const loadAllData = async () => {
        const structureName = await fetchTransplantData(); 
        if (structureName) {
          await fetchStructureData(structureName);
        }
      };
  
      loadAllData();
    }

    if (!window.molstar) {
      console.error("Mol* viewer not found");
      return;
    }
    
    window.molstar.Viewer.create(viewerContainerRef.current!, {
      layoutIsExpanded: false,
      layoutShowControls: false,
      disabledExtensions: ['volseg'],
    }).then(viewer => {
      viewerRef.current = viewer;
    });

  }, [id, best]);

  const handleLoadData = () => {
    if (!viewerRef.current) {
      console.error("Viewer not yet created!");
      return;
    }
  
    const builder = window.molstar.PluginExtensions.mvs.MVSData.createBuilder();
    const structure = builder.download({ url: structure_url }).parse({ format: "mmcif" }).modelStructure({});
    structure.component({ selector: "polymer" }).representation({ type: "cartoon" }).color({ color: "green" });
    const mvsData = builder.getState();
    const mvsj = window.molstar.PluginExtensions.mvs.MVSData.toMVSJ(mvsData);
  
    viewerRef.current.loadMvsData(mvsj, "mvsj", { replaceExisting: true });
  };

  const handleLoadData2 = () => {
    if (!viewerRef.current) {
      console.error("Viewer not yet created!");
      return;
    }
  
    const builder = window.molstar.PluginExtensions.mvs.MVSData.createBuilder();
    const structure = builder.download({ url: structure_url }).parse({ format: "mmcif" }).modelStructure({});
    structure.component({ selector: "polymer" }).representation({ type: "cartoon" }).color({ color: "green" });
    const mvsData = builder.getState();
    const mvsj = window.molstar.PluginExtensions.mvs.MVSData.toMVSJ(mvsData);
  
    viewerRef.current.loadMvsData(mvsj, "mvsj", { replaceExisting: true });
  };

  return (
    <Grid>
      <Grid.Col span={6}>
        {structure && <StructureInfo structure={structure} best={best} setBest={setBest} />}
        <ScrollArea h={"80vh"}>
          {loading ? <Loader /> : 
            <TransplantTable 
              data={transplants} 
              viewerInstanceRef={viewerRef}
              // cognateLigandsData={cognateLigands} 
              // updateCognateLigands={updateCognateLigands} 
            />
          }
          {/* <CognateLigandsTable cognateLigandsData={cognateLigands} /> */}
        </ScrollArea>
      </Grid.Col>
      <Grid.Col span={6}>
        <div ref={viewerContainerRef} style={{ 
              width: "100%", 
              height: "100%", 
              minWidth: "40vh", 
              minHeight: "80vh", 
              position: "relative" }} />
        {/* <MolStarWrapper fileData={structure_url} /> */}
      </Grid.Col>
      <button onClick={handleLoadData}>Load Structure</button>
    </Grid>
  );
}

export default StructurePage;
        //structure?.url