import { useMemo, useState } from 'react';
import {
  MantineReactTable,
  MRT_RowSelectionState,
  useMantineReactTable,
  type MRT_ColumnDef
} from 'mantine-react-table';
import { MantineTheme } from '@mantine/core';
import type { TransplantRecord, CognateLigand } from '../pages/StructurePage';

type FlattenedRecord = CognateLigand & {
  transplant_id: number;
  transplant_ligand_id: number;
  transplant_structure_name: string;
  transplant_ligand_chain: string;
  transplant_tcs: number;
  transplant_name: string;
};

interface TransplantTableProps {
  data: TransplantRecord[]; // Accept data as a prop
  viewerInstanceRef: React.RefObject<any>;
  // cognateLigandsData: CognateLigand[]; // New prop for cognate ligands data
  // updateCognateLigands: (ligands: CognateLigand[]) => void; // Callback to update cognate ligands in the parent
}
// updateCognateLigands
const TransplantTable = ({ data, viewerInstanceRef,  }: TransplantTableProps) => { 
  // table state
  const [currentlyClicked, setCurrentlyClicked] = useState<string | null>(null);
  const [rowSelection, setRowSelection] = useState<MRT_RowSelectionState>({});

  const flattened: FlattenedRecord[] = useMemo(() =>
    data.flatMap((record: TransplantRecord) =>
      record.cognate_ligands.map(cognate => ({
        ...cognate,
        transplant_id: record.transplant.transplant_id,
        transplant_ligand_id: record.transplant.ligand_id,
        transplant_structure_name: record.transplant.structure_name,
        transplant_name: record.transplant.name,
        transplant_ligand_chain: record.transplant.ligand_chain,
        transplant_tcs: record.transplant.tcs
      }))
    ), [data]);

  const columns = useMemo<MRT_ColumnDef<FlattenedRecord>[]>(
    () => [
      {
        accessorKey: 'name',
        header: 'Ligand',
      },
      {
        accessorKey: 'similarity',
        header: 'Similarity'
      },
      {
        accessorKey: 'transplant_tcs',
        header: 'tcs',
      },
      {
        accessorKey: 'transplant_name',
        header: 'Transplanted Ligand'
      },
      {
        accessorKey: 'transplant_ligand_chain',
        header: 'Chain',
      }
    ],
    [],
  );

  const table = useMantineReactTable<FlattenedRecord>({
    columns,
    data : flattened,
    enableMultiRowSelection: false,
    enableFilters: false,
    enableStickyHeader: true,
    enableTopToolbar: false,
    initialState: { showColumnFilters: true },
    mantineTableBodyRowProps: ({ row }) => ({
      onClick: () => {
        const selectedId = row.original.transplant_ligand_chain;
        const alreadySelected = currentlyClicked === selectedId;
        console.log('Selected ID:', selectedId);
        // Update the cognate ligands in the parent component
        if (alreadySelected) {
          viewerInstanceRef.current?.visual.reset({ camera: true });
          setCurrentlyClicked(null);
          setRowSelection({});
        } else {
          viewerInstanceRef.current?.visual.select({
            data: [{ auth_asym_id: selectedId, color: { r: 255, g: 255, b: 0 }, alpha: 0.5, focus: true }],
          });
          setCurrentlyClicked(selectedId);
          setRowSelection({ [row.id]: true });
          console.log('Selected row:', row.original);
          //updateCognateLigands(row.original.cognate_ligands);
        }
      },
      selected: !!rowSelection[row.id], // tells MRT the row is "selected"
      sx: (theme: MantineTheme) => ({
        cursor: 'pointer',
        backgroundColor: rowSelection[row.id]
          ? theme.colors.yellow[1] // or any selected color you want
          : undefined,
      }),
    }),
    state: {
      rowSelection,
    }
  });

  return (
    <div style={{ maxHeight: '80vh', maxWidth: '100%', overflow: 'auto' }}>
      <MantineReactTable table={table} />
    </div>
  );
};

export default TransplantTable;
