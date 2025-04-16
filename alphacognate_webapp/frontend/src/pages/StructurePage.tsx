import { useParams } from "react-router-dom";
import { useEffect, useState, useRef } from "react";
import MolStarViewer from "../components/MolStarViewer";
import TransplantTable from "../components/TransplantTable";
import { Grid, Loader, ScrollArea } from "@mantine/core";
import CognateLigandsTable from "../components/CognateLigandsTable";

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
  const viewerInstanceRef = useRef<any>(null);

  // const [cognateLigands, setCognateLigands] = useState<CognateLigand[]>([]); // State for cognate ligands

  // // Effect to log whenever cognateLigands changes
  // useEffect(() => {
  //   if (cognateLigands.length > 0) {
  //     console.log('Cognate ligands updated:', cognateLigands);
  //   } else {
  //     console.log('Cognate ligands is empty.');
  //   }
  // }, [cognateLigands]); // Dependency array ensures this runs only when cognateLigands changes


  // Callback to update cognate ligands
  // const updateCognateLigands = (ligands: CognateLigand[]) => {
  //   setCognateLigands(ligands); // Update state with new cognate ligands
  // };

  useEffect(() => {
    if (id) {
      const fetchTransplantData = async () => {
        const url = new URL(
          '/transplants',
          process.env.NODE_ENV === 'production'
            ? 'http://localhost:8000/' //replace with prod url when it exists
            : 'http://localhost:8000/',
        );
        url.searchParams.append('structure_name', id);
        try {
          const response = await fetch(url.href);
          const json = (await response.json()) as TransplantApiResponse;
          if (!response.ok) {
            throw new Error('Failed to fetch data');
          }
          setTransplants(json.data); // Update state with fetched data
          // setStructure(json.structure_data);
          let structure_url = new URL(
            json.structure_data?.url,
            process.env.NODE_ENV === 'production'
              ? 'http://localhost:8080/' //replace with prod url when it exists
              : 'http://localhost:8080/',)
          setStructure(json.structure_data);
          setStructureUrl(structure_url.href);
          setLoading(false); // Set loading to false after data is fetched
        } catch (error) {
          console.error(error);
          setLoading(false); // Set loading to false even if an error occurs
        }
      };
      fetchTransplantData();
    }
  }, [id]);

  return (
    <Grid>
      <Grid.Col span={6}>
        <ScrollArea h={"80vh"}>
          {loading ? <Loader /> : 
            <TransplantTable 
              data={transplants} 
              viewerInstanceRef={viewerInstanceRef} 
              // cognateLigandsData={cognateLigands} 
              // updateCognateLigands={updateCognateLigands} 
            />
          }
          {/* <CognateLigandsTable cognateLigandsData={cognateLigands} /> */}
        </ScrollArea>
      </Grid.Col>
      <Grid.Col span={6}>

        <MolStarViewer url={ structure_url || ""} viewerInstanceRef={viewerInstanceRef} /> 
      </Grid.Col>
    </Grid>
  );
}

export default StructurePage;
        //structure?.url