import { useMemo, useState } from 'react';
import {
  MantineReactTable,
  useMantineReactTable,
  type MRT_ColumnDef,
  type MRT_ColumnFiltersState,
  type MRT_PaginationState,
  type MRT_SortingState,
} from 'mantine-react-table';
import type { Transplant } from '../pages/StructurePage';

interface TransplantTableProps {
  data: Transplant[]; // Accept data as a prop
  viewerInstanceRef: React.MutableRefObject<any>;
}

const TransplantTable = ({ data, viewerInstanceRef }: TransplantTableProps) => {
  // table state
  const [columnFilters, setColumnFilters] = useState<MRT_ColumnFiltersState>([]);
  const [globalFilter, setGlobalFilter] = useState('');
  const [currentlyClicked, setCurrentlyClicked] = useState<string | null>(null);
  const [sorting, setSorting] = useState<MRT_SortingState>([]);
  const [pagination, setPagination] = useState<MRT_PaginationState>({
    pageIndex: 0,
    pageSize: 10,
  });

  const columns = useMemo<MRT_ColumnDef<Transplant>[]>(
    () => [
      {
        accessorKey: 'id',
        header: 'Transplant ID',
      },
      {
        accessorKey: 'structure_name',
        header: 'Structure Name',
      },
      {
        accessorKey: 'ligand',
        header: 'Ligand',
      },
      {
        accessorKey: 'tcs',
        header: 'tcs',
      },
      {
        accessorKey: 'struct_asym_id',
        header: 'Chain',
      }
    ],
    [],
  );

  const table = useMantineReactTable({
    columns,
    data,
    enableRowSelection: true,
    enableMultiRowSelection: false,
    initialState: { showColumnFilters: true },
    mantineTableBodyRowProps: ({ row }) => {
      return {
        onClick: (event) => {
          const selectedId = row.original.struct_asym_id;
          if (currentlyClicked === selectedId) {
            // Same row clicked again → Perform a different action
            console.log("Clicked again, performing alternate action:", selectedId);
    
            if (viewerInstanceRef.current) {
              viewerInstanceRef.current.visual.reset({camera: true }); // Example: Clear selection
            }
            
            setCurrentlyClicked(null); // Optionally reset the selection

          } else {
            // Different row clicked → Normal behavior
            console.log("New selection:", selectedId);
    
            if (viewerInstanceRef.current) {
              viewerInstanceRef.current.visual.select({
                 data: [{ auth_asym_id: selectedId, color: { r: 255, g: 255, b: 0 }, focus: true }],
               });
            }
    
            setCurrentlyClicked(selectedId);
          };
          row.getToggleSelectedHandler()(event);
        },
      };
    },
    
    state: {
      columnFilters,
      globalFilter,
      pagination,
      sorting,
    },
  });

  return (
    <div style={{ maxHeight: '80vh', maxWidth: '100%', overflow: 'auto' }}>
      <MantineReactTable table={table} />
    </div>
  );
};

export default TransplantTable;
